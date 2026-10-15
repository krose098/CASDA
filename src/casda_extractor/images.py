"""CASDA image download and cutout generation."""

from __future__ import annotations

import os
import re
import warnings

import astropy.units as u
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.time import Time
from astropy.visualization.wcsaxes import add_scalebar
from astropy.wcs import WCS
from astropy.wcs.wcs import FITSFixedWarning
from astroquery.casda import Casda
from astroquery.utils.tap.core import TapPlus

from .download import parallel_download, stage_in_chunks
from .output import ensure_dir
from .utils import apply_proper_motion

warnings.filterwarnings("ignore", category=FITSFixedWarning)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _simplify_fits(hdu) -> tuple:
    """
    Collapse 3-D or 4-D ASKAP FITS data to 2-D and strip extra WCS axes.

    Returns
    -------
    (data_2d, header_2d)
    """
    shape = hdu.data.shape
    if len(shape) == 4:
        data = hdu.data[0, 0, :, :]
    elif len(shape) == 3:
        data = hdu.data[0, :, :]
    else:
        data = hdu.data.copy()

    header = hdu.header.copy()
    remove = ["NAXIS3", "NAXIS4"]
    for axis in (3, 4):
        for kw in ["CRPIX", "CTYPE", "CRVAL", "CDELT", "CUNIT"]:
            remove.append(f"{kw}{axis}")
        for i in range(1, 5):
            for j in range(1, 5):
                if i > 2 or j > 2:
                    remove += [f"PC{i}_{j}", f"CD{i}_{j}"]
    for k in remove:
        if k in header:
            del header[k]
    header["NAXIS"] = 2
    if "WCSAXES" in header:
        header["WCSAXES"] = 2
    return data, header


def _rebuild_cutout_metadata(source_name: str, cutout_dir: str) -> pd.DataFrame:
    """Reconstruct cutout metadata by scanning an existing cutout directory."""
    clean = source_name.replace(" ", "").replace("(", "").replace(")", "")
    import glob

    files = glob.glob(os.path.join(cutout_dir, f"{clean}_*.fits"))
    records = []
    for fp in files:
        try:
            with fits.open(fp) as hdul:
                hdr = hdul[0].header
                obs_date = hdr.get("OBSDATE", hdr.get("DATE-OBS", ""))
                mjd = hdr.get("MJD-OBS", 0)
                if mjd == 0 and obs_date:
                    try:
                        mjd = Time(obs_date).mjd
                    except Exception:  # noqa: BLE001
                        mjd = 0
                records.append(
                    {
                        "source_name": source_name,
                        "cutout_filename": fp,
                        "image_filename": os.path.basename(fp),
                        "obs_date": obs_date,
                        "obs_mjd": mjd,
                        "integration_time": hdr.get("INT_TIME", 0),
                    }
                )
        except Exception:  # noqa: BLE001
            pass
    return pd.DataFrame(records)


def _filter_valid_cutouts(df: pd.DataFrame) -> pd.DataFrame:
    """Drop cutout records whose files are all-NaN or unreadable."""
    valid = []
    for _, row in df.iterrows():
        try:
            with fits.open(row["cutout_filename"]) as hdul:
                if np.any(~np.isnan(hdul[0].data)):
                    valid.append(row)
        except Exception:  # noqa: BLE001
            pass
    return pd.DataFrame(valid)


# ---------------------------------------------------------------------------
# Query
# ---------------------------------------------------------------------------


def query_images(
    ra: float,
    dec: float,
    search_radius: float = 10.0,
    image_type: str = "cont.restored.t0",
    stokes: str = "I",
) -> pd.DataFrame:
    """
    Query CASDA TAP for public ASKAP continuum images near ``(ra, dec)``.

    Parameters
    ----------
    ra, dec : float
        ICRS coordinates in degrees.
    search_radius : float
        Search cone in degrees.
    image_type : str
        CASDA ``dataproduct_subtype`` value.
    stokes : str
        Stokes parameter filter.

    Returns
    -------
    pd.DataFrame
        Filtered image metadata.
    """
    tap = TapPlus(url="https://casda.csiro.au/casda_vo_tools/tap")
    print(f"Querying CASDA for Stokes {stokes} images (type: {image_type})…")
    job = tap.launch_job_async(
        f"SELECT TOP 50000 * FROM ivoa.obscore "
        f"WHERE pol_states = '/{stokes}/' AND dataproduct_subtype = '{image_type}'"
    )
    r = job.get_results()
    data = r[(r["quality_level"] == "GOOD") | (r["quality_level"] == "UNCERTAIN")]
    public_data = Casda.filter_out_unreleased(data).to_pandas()

    # Keep only ASKAP / Beta observations
    askap_mask = (
        public_data["obs_id"].str.contains("ASKAP", na=False)
        | public_data["obs_id"].str.contains("BETA", na=False)
        | public_data["obs_id"].str.contains("Beta", na=False)
    )
    public_data = public_data[askap_mask]
    print(f"Total public ASKAP images: {len(public_data)}")

    # Spatial filter
    src = SkyCoord(ra * u.deg, dec * u.deg)
    img_coords = SkyCoord(
        np.array(public_data["s_ra"]) * u.deg,
        np.array(public_data["s_dec"]) * u.deg,
    )
    matches = np.where(src.separation(img_coords) < search_radius * u.deg)[0]
    matched = public_data.iloc[matches].copy()
    print(f"Images covering source: {len(matched)}")

    # Suffix filter
    suffix = ".restored.conv.fits"
    matched = matched[matched["filename"].str.endswith(suffix)].copy()
    print(f'After suffix filter ("{suffix}"): {len(matched)} images')
    return matched


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------


def download_images(
    image_data: pd.DataFrame,
    image_dir: str,
    casda_obj,
    parallel: bool = True,
    max_workers: int = 8,
    chunk_size: int = 50,
) -> list[str]:
    """
    Download CASDA images, skipping already-downloaded files.

    Parameters
    ----------
    image_data : pd.DataFrame
        Metadata from ``query_images()``.
    image_dir : str
        Directory to save images.
    casda_obj : Casda
        Authenticated Casda instance.
    parallel : bool
        Use parallel downloads.
    max_workers : int
        Worker threads for parallel mode.
    chunk_size : int
        Staging chunk size.

    Returns
    -------
    list[str]
        Paths of downloaded files.
    """
    ensure_dir(image_dir)
    if not image_dir.endswith("/"):
        image_dir += "/"

    to_download = [f for f in image_data["filename"] if not os.path.isfile(os.path.join(image_dir, f))]
    if not to_download:
        print("All images already downloaded.")
        return []

    print(f"Images to download: {len(to_download)}")
    url_list = stage_in_chunks(to_download, image_data, casda_obj, chunk_size=chunk_size)
    print(f"Staged {len(url_list)} URLs.")

    if parallel and len(url_list) > 1:
        return parallel_download(url_list, image_dir, casda_obj, max_workers=max_workers)
    return casda_obj.download_files(url_list, savedir=image_dir)


# ---------------------------------------------------------------------------
# Cutout generation
# ---------------------------------------------------------------------------


def make_cutouts(
    sources_df: pd.DataFrame,
    image_data: pd.DataFrame,
    image_dir: str,
    cutout_dir: str,
    cutout_size: u.Quantity = 4.0 * u.arcmin,
    ra_col: str = "ra",
    dec_col: str = "dec",
    pmra_col: str = "pmra",
    pmdec_col: str = "pmdec",
    epoch_col: str = "epoch",
    name_col: str = "name",
    apply_pm: bool = True,
) -> pd.DataFrame:
    """
    Generate per-source cutouts from CASDA full-field images.

    Cutouts that are entirely NaN are silently skipped.

    Parameters
    ----------
    sources_df : pd.DataFrame
        Source table with at minimum ``name``, ``ra``, ``dec`` columns.
    image_data : pd.DataFrame
        Metadata from ``query_images()``.
    image_dir : str
        Directory containing the full-field FITS images.
    cutout_dir : str
        Directory to write cutout FITS files.
    cutout_size : astropy.units.Quantity
        Angular size of the cutout square.
    ra_col, dec_col : str
        Column names for RA and Dec in ``sources_df``.
    pmra_col, pmdec_col : str
        Proper-motion column names (mas/yr).
    epoch_col : str
        Reference epoch column name (e.g. "J2000").
    name_col : str
        Source name column.
    apply_pm : bool
        Whether to propagate proper motion to the image epoch.

    Returns
    -------
    pd.DataFrame
        Records of successfully created cutouts.
    """
    ensure_dir(cutout_dir)
    if not image_dir.endswith("/"):
        image_dir += "/"
    if not cutout_dir.endswith("/"):
        cutout_dir += "/"

    records = []

    for _, img_row in image_data.iterrows():
        filename = img_row["filename"]
        filepath = image_dir + filename
        if not os.path.isfile(filepath):
            continue

        t_min_val = img_row["t_min"]
        if not np.isfinite(t_min_val):
            print(f"  Skipping {filename}: invalid t_min.")
            continue

        img_epoch = Time(t_min_val, format="mjd")
        img_coord = SkyCoord(img_row["s_ra"] * u.deg, img_row["s_dec"] * u.deg)

        # Only process sources within ~6° of image centre
        nearby = [
            src
            for _, src in sources_df.iterrows()
            if SkyCoord(src[ra_col] * u.deg, src[dec_col] * u.deg).separation(img_coord) < 6 * u.deg
        ]
        if not nearby:
            continue

        print(f"\nProcessing: {filename}")
        hdul = None
        try:
            hdul = fits.open(filepath)
            hdu = hdul[0]
            original_header = hdu.header
            data, header = _simplify_fits(hdu)
            wcs = WCS(header)
            pix_scale = abs(header["CDELT1"]) * u.deg
            cutout_pixels = int((cutout_size / pix_scale).decompose().value)

            for src in nearby:
                src_name = src[name_col]
                clean = src_name.replace(" ", "").replace("(", "").replace(")", "")
                cutout_path = f"{cutout_dir}{clean}_{filename}"
                if os.path.isfile(cutout_path):
                    print(f"  Skipping existing: {clean}")
                    continue

                ra_val, dec_val = float(src[ra_col]), float(src[dec_col])
                src_coord = SkyCoord(ra_val * u.deg, dec_val * u.deg)

                if apply_pm and pmra_col in src and pmdec_col in src:
                    pmra_v = float(np.nan_to_num(src[pmra_col], nan=0.0))
                    pmdec_v = float(np.nan_to_num(src[pmdec_col], nan=0.0))
                    src_coord = apply_proper_motion(
                        src_coord, pmra_v, pmdec_v, str(src[epoch_col]), t_min_val
                    )

                try:
                    cutout = Cutout2D(data, src_coord, size=cutout_pixels, wcs=wcs)
                    if not np.any(~np.isnan(cutout.data)):
                        print(f"  Skipping all-NaN cutout: {src_name}")
                        continue

                    chdu = fits.PrimaryHDU(data=cutout.data)
                    chdu.header.update(cutout.wcs.to_header())
                    for k in ["BMAJ", "BMIN", "BPA", "BUNIT", "TELESCOP", "INSTRUME"]:
                        if k in original_header:
                            chdu.header[k] = original_header[k]
                    chdu.header["OBJECT"] = src_name
                    chdu.header["OBSDATE"] = img_epoch.isot
                    chdu.header["MJD-OBS"] = img_row["t_min"]
                    chdu.header["INT_TIME"] = img_row["t_exptime"]
                    chdu.writeto(cutout_path, overwrite=True)

                    records.append(
                        {
                            "source_name": src_name,
                            "image_filename": filename,
                            "cutout_filename": cutout_path,
                            "obs_date": img_epoch.isot,
                            "obs_mjd": img_row["t_min"],
                            "integration_time": img_row["t_exptime"],
                        }
                    )
                    print(f"  Created: {clean}")
                except Exception as exc:  # noqa: BLE001
                    print(f"  Failed cutout {src_name}: {exc}")

        except Exception as exc:  # noqa: BLE001
            print(f"Error processing {filename}: {exc}")
        finally:
            if hdul:
                hdul.close()

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------


def plot_cutout_grid(
    source_name: str,
    cutout_info: pd.DataFrame,
    cutout_dir: str,
    save_path: str | None = None,
    target_ra: float | None = None,
    target_dec: float | None = None,
) -> str:
    """
    Save a grid of cutout images with WCS axes, beam ellipses, and optional crosshair.

    Parameters
    ----------
    source_name : str
        Source label (used to filter ``cutout_info`` and in the title).
    cutout_info : pd.DataFrame
        Valid cutout records from ``make_cutouts()`` or ``_rebuild_cutout_metadata()``.
    cutout_dir : str
        Directory containing the cutout FITS files.
    save_path : str, optional
        Directory for the output PNG (defaults to ``cutout_dir``).
    target_ra, target_dec : float, optional
        If given, a crosshair is drawn at this sky position.

    Returns
    -------
    str
        Path to the saved PNG.
    """
    if not cutout_dir.endswith("/"):
        cutout_dir += "/"
    if save_path is None:
        save_path = cutout_dir
    elif not save_path.endswith("/"):
        save_path += "/"

    src_cuts = cutout_info[cutout_info["source_name"] == source_name].sort_values("obs_mjd")
    n = len(src_cuts)
    if n == 0:
        print(f"No valid cutouts for {source_name}.")
        return ""

    n_cols, n_rows = 3, (n + 2) // 3
    fig = plt.figure(figsize=(12, 3.8 * n_rows))

    for idx, (_, row) in enumerate(src_cuts.iterrows()):
        try:
            with fits.open(row["cutout_filename"]) as hdul:
                hdr = hdul[0].header
                data = hdul[0].data
                wcs = WCS(hdr)
                bmaj = hdr.get("BMAJ", 0.0)
                bmin = hdr.get("BMIN", 0.0)
                bpa = hdr.get("BPA", 0.0)
                pix_scale = abs(hdr.get("CDELT2", hdr.get("CD2_2", 0.0)))

            ax = fig.add_subplot(n_rows, n_cols, idx + 1, projection=wcs)
            vmin, vmax = np.nanpercentile(data[~np.isnan(data)], [1, 99])
            im = ax.imshow(data, cmap="magma", origin="lower", vmin=vmin, vmax=vmax)
            cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            cb.set_label(hdr.get("BUNIT", "Jy/beam"), fontsize=8)
            cb.ax.tick_params(labelsize=7)

            ax.coords[0].set_axislabel("RA (J2000)", fontsize=8)
            ax.coords[1].set_axislabel("Dec (J2000)", fontsize=8)
            ax.coords[0].set_ticklabel(size=7)
            ax.coords[1].set_ticklabel(size=7)
            ax.coords[0].set_major_formatter("hh:mm:ss.s")
            ax.coords[1].set_major_formatter("dd:mm:ss")

            if target_ra is not None and target_dec is not None:
                tc = SkyCoord(target_ra * u.deg, target_dec * u.deg)
                ax.plot_coord(tc, "w+", markersize=12, markeredgewidth=2, markeredgecolor="gray", zorder=15)

            if bmaj > 0 and pix_scale > 0:
                maj_pix = bmaj / pix_scale
                min_pix = bmin / pix_scale
                bx, by = data.shape[1] * 0.1, data.shape[0] * 0.1
                ax.add_patch(
                    patches.Ellipse(
                        (bx, by),
                        width=min_pix,
                        height=maj_pix,
                        angle=bpa,
                        edgecolor="cyan",
                        facecolor="none",
                        hatch="///",
                        linewidth=1.5,
                        zorder=10,
                    )
                )
                ax.set_xlim(-0.5, data.shape[1] - 0.5)
                ax.set_ylim(-0.5, data.shape[0] - 0.5)
                cx, cy = data.shape[1] // 2, data.shape[0] // 2
                ax.add_patch(
                    patches.Rectangle(
                        (cx - maj_pix / 2, cy - maj_pix / 2),
                        maj_pix,
                        maj_pix,
                        linewidth=1,
                        edgecolor="white",
                        facecolor="none",
                    )
                )

            obs_date = str(row["obs_date"]).split("T")[0]
            int_time = round(row["integration_time"] / 60, 1)
            beam_str = f"{bmaj * 3600:.2f} × {bmin * 3600:.2f} ''"
            sbid_match = re.search(r"(SB\d+)", row["image_filename"])
            sbid_str = sbid_match.group(1) if sbid_match else ""
            ax.set_title(f"{obs_date} | {sbid_str}\n{int_time} min\nbeam: {beam_str}", fontsize=9)
            add_scalebar(ax, 5 * u.arcsec, label="5''", color="white", corner="top right")

        except Exception as exc:  # noqa: BLE001
            ax = fig.add_subplot(n_rows, n_cols, idx + 1)
            ax.text(0.5, 0.5, "Error", ha="center", va="center", transform=ax.transAxes)
            ax.set_xticks([])
            ax.set_yticks([])
            print(f"Plot error {row['image_filename']}: {exc}")

    fig.suptitle(source_name, fontsize=14, y=0.99)
    fig.tight_layout()

    clean = source_name.replace(" ", "").replace("(", "").replace(")", "")
    out_path = os.path.join(save_path, f"{clean}_cutouts.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Cutout grid saved: {out_path}")
    return out_path


# ---------------------------------------------------------------------------
# High-level entry point
# ---------------------------------------------------------------------------


def run_images(
    name: str,
    ra: float,
    dec: float,
    output_root: str,
    casda_obj,
    pmra: float = 0.0,
    pmdec: float = 0.0,
    epoch: str = "J2000",
    stokes: str = "I",
    search_radius: float = 10.0,
    cutout_size: float = 4.0,
    parallel: bool = True,
    max_workers: int = 8,
    plot_only: bool = False,
) -> pd.DataFrame:
    """
    Full images + cutouts pipeline for a single source.

    Parameters
    ----------
    name : str
        Source label.
    ra, dec : float
        ICRS coordinates in degrees.
    output_root : str
        Root output directory.  Files written under ``output_root/{name}/``.
    casda_obj : Casda
        Authenticated Casda instance (ignored when ``plot_only=True``).
    pmra, pmdec : float
        Proper motion in mas/yr (set both to 0 to disable).
    epoch : str
        Reference epoch string for proper motion (e.g. "J2000").
    stokes : str
        Stokes parameter.
    search_radius : float
        Image search cone in degrees.
    cutout_size : float
        Cutout size in arcminutes.
    parallel : bool
        Parallel image downloads.
    max_workers : int
        Download thread count.
    plot_only : bool
        Skip query/download and only (re-)plot existing cutouts.

    Returns
    -------
    pd.DataFrame
        Cutout metadata.
    """
    src_dir = os.path.join(output_root, name)
    image_dir = os.path.join(src_dir, "images")
    cutout_dir = os.path.join(src_dir, "cutouts")
    ensure_dir(image_dir)
    ensure_dir(cutout_dir)

    if plot_only:
        print(f"[{name}] Plot-only mode — scanning existing cutouts…")
        cutout_info = _rebuild_cutout_metadata(name, cutout_dir)
    else:
        sources_df = pd.DataFrame(
            {"name": [name], "ra": [ra], "dec": [dec], "pmra": [pmra], "pmdec": [pmdec], "epoch": [epoch]}
        )
        image_data = query_images(ra, dec, search_radius=search_radius, stokes=stokes)
        if image_data.empty:
            print(f"No images found for {name}.")
            return pd.DataFrame()

        download_images(image_data, image_dir, casda_obj, parallel=parallel, max_workers=max_workers)
        cutout_info = make_cutouts(
            sources_df,
            image_data,
            image_dir,
            cutout_dir,
            cutout_size=cutout_size * u.arcmin,
            apply_pm=(pmra != 0.0 or pmdec != 0.0),
        )

    if cutout_info is None or cutout_info.empty:
        print(f"No cutouts for {name}.")
        return pd.DataFrame()

    cutout_info = _filter_valid_cutouts(cutout_info)
    if cutout_info.empty:
        print(f"No valid cutouts for {name} after filtering.")
        return pd.DataFrame()

    plot_cutout_grid(name, cutout_info, cutout_dir, target_ra=ra, target_dec=dec)
    print(f"[{name}] Done — {len(cutout_info)} valid cutouts.")
    return cutout_info
