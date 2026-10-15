"""Detection lightcurve extraction from CASDA catalogue XML files."""

from __future__ import annotations

import os
import re
from datetime import datetime

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io.votable import parse
from astroquery.casda import Casda
from astroquery.utils.tap.core import TapPlus

from .download import parallel_download, stage_in_chunks
from .output import ensure_dir
from .utils import apply_proper_motion

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def _convert_xml_to_df(xml_path: str) -> pd.DataFrame:
    """Parse a CASDA VOTable XML file into a pandas DataFrame."""
    votable = parse(xml_path)
    table = votable.get_first_table()
    return table.to_table(use_names_over_ids=True).to_pandas()


def _extract_sbid(filename: str) -> str | None:
    """Return the numeric SBID string from a filename, or None."""
    match = re.search(r"SB(\d+)", filename)
    return match.group(1) if match else None


def _error_quad(flux: float, flux_err: float, rms: float, scale: float = 0.06) -> float:
    """Quadrature error combining fit uncertainty, image RMS, and flux scale error."""
    return np.sqrt(flux_err**2 + rms**2 + (scale * flux) ** 2)


# ---------------------------------------------------------------------------
# Download catalogues
# ---------------------------------------------------------------------------


def download_catalogues(
    name: str,
    ra: float,
    dec: float,
    outdir: str,
    casda_obj,
    search_radius: float = 10.0,
    max_workers: int = 4,
    chunk_size: int = 50,
) -> pd.DataFrame:
    """
    Query CASDA for continuum catalogues near ``(ra, dec)`` and download them.

    Parameters
    ----------
    name : str
        Source label (used for the save subdirectory).
    ra, dec : float
        ICRS coordinates in degrees.
    outdir : str
        Root output directory for this source (catalogues saved to ``outdir/catalogues/``).
    casda_obj : Casda
        Authenticated Casda instance.
    search_radius : float
        Search cone in degrees.
    max_workers : int
        Parallel download threads.
    chunk_size : int
        Staging chunk size.

    Returns
    -------
    pd.DataFrame
        Public catalogue metadata from CASDA TAP.
    """
    cat_dir = os.path.join(outdir, "catalogues")
    ensure_dir(cat_dir)

    tap = TapPlus(url="https://casda.csiro.au/casda_vo_tools/tap")
    print("Querying CASDA TAP for continuum catalogues…")
    job = tap.launch_job_async(
        "SELECT TOP 50000 * FROM ivoa.obscore " "WHERE dataproduct_subtype = 'catalogue.continuum.component'"
    )
    r = job.get_results()
    data = r[(r["quality_level"] == "GOOD") | (r["quality_level"] == "UNCERTAIN")]
    public_data = Casda.filter_out_unreleased(data)
    pubdat = public_data.to_pandas()

    print(f"Total public catalogues in archive: {len(pubdat)}")
    print(f"Source: RA={ra:.4f}°, Dec={dec:.4f}°  |  Search radius: {search_radius}°")

    src_coord = SkyCoord(ra * u.deg, dec * u.deg)
    img_coords = SkyCoord(
        np.array(public_data["s_ra"]) * u.deg,
        np.array(public_data["s_dec"]) * u.deg,
    )
    seps = src_coord.separation(img_coords)
    matches_idx = np.where(seps < search_radius * u.deg)[0]
    matching_files = [f for f in np.array(pubdat.iloc[matches_idx]["filename"]) if "checksum" not in f]
    print(f"Matching catalogue files: {len(matching_files)}")

    if not matching_files:
        return pubdat

    # Stage and download
    print("Staging files…")
    url_list = stage_in_chunks(matching_files, pubdat, casda_obj, chunk_size=chunk_size)
    print(f"Staged {len(url_list)} URLs. Downloading…")

    if len(url_list) > 1:
        parallel_download(url_list, cat_dir, casda_obj, max_workers=max_workers)
    else:
        casda_obj.download_files(url_list, savedir=cat_dir)

    print("Catalogue download complete.")
    return pubdat


# ---------------------------------------------------------------------------
# Filter lightcurve matches
# ---------------------------------------------------------------------------


def build_lightcurve(
    name: str,
    ra: float,
    dec: float,
    outdir: str,
    pubdat: pd.DataFrame,
    sep: float = 5.0,
    scale: float = 0.06,
    pmra: float = 0.0,
    pmdec: float = 0.0,
    epoch: str = "J2000",
) -> tuple[pd.DataFrame, list]:
    """
    Parse downloaded catalogues, cross-match to source, and build a lightcurve DataFrame.

    Parameters
    ----------
    name : str
        Source label.
    ra, dec : float
        ICRS coordinates in degrees.
    outdir : str
        Source output directory (catalogues expected in ``outdir/catalogues/``).
    pubdat : pd.DataFrame
        Metadata from ``download_catalogues()`` (used to attach obs times).
    sep : float
        Cross-match radius in arcseconds.
    scale : float
        Fractional flux-scale uncertainty for quadrature error.
    pmra, pmdec : float
        Proper motion in mas/yr (RA × cos δ and Dec).  Set both to 0 to disable.
    epoch : str
        Reference epoch for the catalogue position (e.g. ``"J2000"`` or ``"2016.0"``).

    Returns
    -------
    (df_lc, no_match_list)
        ``df_lc`` sorted lightcurve DataFrame; ``no_match_list`` epochs with no detection.
    """
    cat_dir = os.path.join(outdir, "catalogues")
    base_target = SkyCoord(ra * u.deg, dec * u.deg)
    max_sep = sep * u.arcsec
    apply_pm = pmra != 0.0 or pmdec != 0.0

    # Build obs-info lookup from metadata
    obs_lookup: dict[str, dict] = {}
    for _, row in pubdat.iterrows():
        obs_id = row["obs_id"]
        if obs_id not in obs_lookup:
            obs_lookup[obs_id] = {
                "t_min": row["t_min"],
                "t_exptime": row["t_exptime"],
                "obs_collection": row["obs_collection"],
            }

    data_rows = []
    no_match_list = []

    for fname in os.listdir(cat_dir):
        if not fname.endswith(".components.xml"):
            continue
        sbid = _extract_sbid(fname)
        df = _convert_xml_to_df(os.path.join(cat_dir, fname))
        if len(df) == 0:
            no_match_list.append(fname)
            continue

        # Resolve observation info early so we can compute PM-corrected position
        obs_info: dict = {}
        for prefix in ["ASKAP-", "Beta-"]:
            key = f"{prefix}{sbid}"
            if key in obs_lookup:
                obs_info = obs_lookup[key]
                break
        if not obs_info:
            obs_info = {"t_min": np.nan, "t_exptime": np.nan, "obs_collection": "Unknown"}

        # Compute per-epoch target position (proper-motion corrected if requested)
        target = base_target
        if apply_pm:
            t_min = obs_info.get("t_min", np.nan)
            target = apply_proper_motion(base_target, pmra, pmdec, epoch, t_min)

        df_coords = SkyCoord(
            ra=df["ra_deg_cont"].values * u.deg,
            dec=df["dec_deg_cont"].values * u.deg,
        )
        seps = target.separation(df_coords)
        mask = seps <= max_sep

        if not mask.any():
            no_match_list.append(fname)
            continue

        idx = np.where(mask)[0][0]
        row = df.iloc[idx]

        data_rows.append(
            {
                "sbid": sbid,
                "component_id": row["component_id"],
                "component_name": row["component_name"],
                "ra_deg_cont": row["ra_deg_cont"],
                "dec_deg_cont": row["dec_deg_cont"],
                "freq": row["freq"],
                "flux_peak": row["flux_peak"],
                "flux_peak_err": row["flux_peak_err"],
                "rms_image": row["rms_image"],
                "obs_start": obs_info["t_min"],
                "obs_length": obs_info["t_exptime"],
                "project": obs_info["obs_collection"],
            }
        )

    if not data_rows:
        return pd.DataFrame(), no_match_list

    df_lc = pd.DataFrame(data_rows)
    df_lc["flux_err_quad"] = _error_quad(df_lc.flux_peak, df_lc.flux_peak_err, df_lc.rms_image, scale)
    df_lc = df_lc.sort_values("obs_start").reset_index(drop=True)

    print(f"Lightcurve built: {len(df_lc)} detections, {len(no_match_list)} non-detections.")
    return df_lc, no_match_list


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------


def plot_lightcurve(
    df: pd.DataFrame,
    name: str,
    outdir: str,
) -> str:
    """
    Save a scatter plot of the detection lightcurve.

    Parameters
    ----------
    df : pd.DataFrame
        Output of ``build_lightcurve()``.
    name : str
        Source label (used in title and filename).
    outdir : str
        Source output directory (plot saved to ``outdir/``).

    Returns
    -------
    str
        Path to the saved PDF.
    """
    if df.empty:
        print("No detections to plot.")
        return ""

    norm = plt.Normalize(df["freq"].min(), df["freq"].max())
    cmap = plt.get_cmap("viridis_r")

    fig, ax = plt.subplots()
    sc = ax.scatter(df["obs_start"], df["flux_peak"], c=df["freq"], cmap=cmap, norm=norm, s=50)
    ax.errorbar(
        df["obs_start"],
        df["flux_peak"],
        yerr=df["flux_err_quad"],
        fmt="none",
        ecolor="gray",
        capsize=2,
        alpha=0.4,
    )
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("Frequency [MHz]")
    ax.set_xlabel("Date [MJD]")
    ax.set_ylabel("Flux Density [mJy/beam]")
    ax.set_title(name)
    fig.tight_layout()

    timestamp = datetime.now().strftime("%Y%m%d_%H%M")
    save_path = os.path.join(outdir, f"{name}_lc_{timestamp}.pdf")
    fig.savefig(save_path, dpi=300)
    plt.close(fig)
    print(f"Lightcurve plot saved: {save_path}")
    return save_path


# ---------------------------------------------------------------------------
# High-level entry point
# ---------------------------------------------------------------------------


def run_lightcurves(
    name: str,
    ra: float,
    dec: float,
    output_root: str,
    casda_obj,
    search_radius: float = 10.0,
    sep: float = 5.0,
    scale: float = 0.06,
    max_workers: int = 4,
    chunk_size: int = 50,
    plot: bool = True,
    pmra: float = 0.0,
    pmdec: float = 0.0,
    epoch: str = "J2000",
) -> pd.DataFrame:
    """
    Full lightcurve pipeline: query → download → cross-match → plot → save CSV.

    Parameters
    ----------
    name : str
        Source label.
    ra, dec : float
        ICRS coordinates in degrees.
    output_root : str
        Root output directory.  Results written to ``output_root/{name}/lightcurves/``.
    casda_obj : Casda
        Authenticated Casda instance.
    search_radius : float
        TAP cone search radius in degrees.
    sep : float
        Cross-match radius in arcseconds.
    scale : float
        Fractional flux-scale uncertainty (default 6 %).
    max_workers : int
        Parallel download threads.
    chunk_size : int
        Staging chunk size.
    plot : bool
        Whether to generate and save the scatter plot.
    pmra, pmdec : float
        Proper motion in mas/yr.  Set both to 0 to disable.
    epoch : str
        Reference epoch for the catalogue position (e.g. ``"J2000"``).

    Returns
    -------
    pd.DataFrame
        Lightcurve DataFrame (empty if no detections).
    """
    lc_dir = os.path.join(output_root, name, "lightcurves")
    ensure_dir(lc_dir)

    pubdat = download_catalogues(
        name,
        ra,
        dec,
        lc_dir,
        casda_obj,
        search_radius=search_radius,
        max_workers=max_workers,
        chunk_size=chunk_size,
    )
    df_lc, _ = build_lightcurve(
        name, ra, dec, lc_dir, pubdat, sep=sep, scale=scale, pmra=pmra, pmdec=pmdec, epoch=epoch
    )

    if not df_lc.empty:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M")
        csv_path = os.path.join(lc_dir, f"{name}_lc_{timestamp}.csv")
        df_lc.to_csv(csv_path, index=False)
        print(f"Lightcurve CSV saved: {csv_path}")
        if plot:
            plot_lightcurve(df_lc, name, lc_dir)

    return df_lc
