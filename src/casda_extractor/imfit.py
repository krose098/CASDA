"""Forced-fit photometry — CASA imfit/imstat when available, astropy fallback otherwise."""

from __future__ import annotations

import csv
import os
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.time import Time

from .output import ensure_dir
from .utils import apply_proper_motion, extract_fits_header, unicoord

try:
    from casatasks import imfit, imstat  # type: ignore[import]

    _CASA_AVAILABLE = True
except ImportError:
    _CASA_AVAILABLE = False


# ---------------------------------------------------------------------------
# CASA helpers (used only when casatasks is available)
# ---------------------------------------------------------------------------


def _coords_to_casa_region(coords: SkyCoord, radius: float, rms_radius: float) -> tuple[str, str]:
    """
    Build CASA region strings (circular fit region + annulus for RMS).

    Parameters
    ----------
    coords : SkyCoord
        Source position.
    radius : float
        Fit circle radius in arcseconds.
    rms_radius : float
        Outer annulus radius in arcseconds.

    Returns
    -------
    (fit_region, rms_annulus_region)
    """
    hmsdms = coords.to_string("hmsdms", sep=":", precision=1, pad=True)
    ra_str, dec_str = hmsdms.split()
    hh, mm, ss = ra_str.split(":")
    ra_casa = f"{hh}h{mm}m{ss}s"
    sign = "+" if dec_str.strip()[0] == "+" else "-"
    dd, dm, ds = dec_str.replace("+", "").replace("-", "").split(":")
    dec_casa = f"{sign}{dd}.{dm}.{ds[:4]}"
    centre = f"{ra_casa},{dec_casa}"
    fit_region = f"circle[[{centre}], {radius}arcsec]"
    rms_annulus = f"annulus[[{centre}], [{2 * radius}arcsec, {rms_radius}arcsec]]"
    return fit_region, rms_annulus


def _fit_image_casa(
    path: str,
    coords: SkyCoord,
    stokes: str = "I",
    radius: float = 5.0,
    rms_radius: float = 30.0,
    scale: float = 0.1,
    units: str = "mJy",
    workdir: str | None = None,
) -> dict:
    """Run CASA imfit + imstat on a single FITS image."""
    unit_scale = 1e6 if units == "uJy" else 1e3
    unit_label = units

    header = extract_fits_header(path)
    result: dict = {
        "Date of Observation (UTC)": header.get("DATE-OBS"),
        "Project": header.get("PROJECT"),
        "Freq (GHz)": header.get("CRVAL3"),
        f"S ({unit_label})": None,
        f"Err ({unit_label})": None,
        f"rms ({unit_label}/beam)": None,
        "SBID": header.get("SBID"),
        "Duration": header.get("DURATION"),
        "file": path,
        "ra_fit_deg": None,
        "ra_err_deg": None,
        "dec_fit_deg": None,
        "dec_err_deg": None,
        "fit_failed": False,
    }

    fit_region, rms_annulus = _coords_to_casa_region(coords, radius, rms_radius)

    ra_str = coords.to_string("hmsdms", sep=":", precision=1, pad=True).split()[0]
    hh, mm = ra_str.split(":")[:2]
    source_label = f"J{hh}{mm}_{stokes}_{radius}arcsec"
    prefix = os.path.join(workdir or ".", source_label)

    try:
        fit = imfit(
            imagename=path,
            region=fit_region,
            chans="",
            stokes="",
            summary=f"{prefix}_summary",
            logfile=f"{prefix}_log",
        )
        fit_failed = not fit.get("converged", [False])[0]
    except Exception:  # noqa: BLE001
        fit = None
        fit_failed = True

    stats = imstat(imagename=path, stokes=stokes, region=rms_annulus)
    rms = stats["rms"] * unit_scale
    result[f"rms ({unit_label}/beam)"] = round(float(rms[0]), 4)
    result["fit_failed"] = fit_failed

    if not fit_failed and fit is not None:
        comp = fit["results"]["component0"]
        peak_flux = comp["peak"]["value"] * unit_scale
        peak_err = comp["peak"]["error"] * unit_scale
        combined_err = (peak_err**2 + rms**2 + (scale * peak_flux) ** 2) ** 0.5

        direction = comp["shape"]["direction"]
        ra_deg = direction["m0"]["value"] * 180 / np.pi % 360
        ra_err_deg = direction["error"]["longitude"]["value"] / 3600
        dec_deg = direction["m1"]["value"] * 180 / np.pi
        dec_err_deg = direction["error"]["latitude"]["value"] / 3600

        result.update(
            {
                f"S ({unit_label})": round(float(peak_flux), 4),
                f"Err ({unit_label})": round(float(combined_err[0]), 4),
                "ra_fit_deg": round(ra_deg, 6),
                "ra_err_deg": round(ra_err_deg, 6),
                "dec_fit_deg": round(dec_deg, 6),
                "dec_err_deg": round(dec_err_deg, 6),
            }
        )

    return result


# ---------------------------------------------------------------------------
# Pure-Python fallback (astropy + numpy)
# ---------------------------------------------------------------------------


def _fit_image_astropy(
    path: str,
    coords: SkyCoord,
    stokes: str = "I",  # noqa: ARG001  (kept for API symmetry with CASA path)
    radius: float = 5.0,
    rms_radius: float = 30.0,
    scale: float = 0.1,
    units: str = "mJy",
    workdir: str | None = None,  # noqa: ARG001
) -> dict:
    """
    Fit a 2-D Gaussian to a FITS image cutout using astropy.modeling.

    Used automatically when ``casatasks`` is not installed (e.g. on macOS).
    Results are equivalent to CASA imfit for point/slightly-resolved sources.

    Parameters
    ----------
    path : str
        Path to the FITS image (Jy/beam).
    coords : SkyCoord
        Source position.
    stokes : str
        Ignored (kept for signature compatibility with the CASA path).
    radius : float
        Fit-circle radius in arcseconds.
    rms_radius : float
        Outer annulus radius (arcsec) for background RMS estimation.
    scale : float
        Fractional flux-scale uncertainty added in quadrature.
    units : str
        Output flux units: ``'mJy'`` or ``'uJy'``.
    workdir : str, optional
        Ignored (kept for signature compatibility with the CASA path).

    Returns
    -------
    dict
        Result row matching the schema produced by the CASA path.
    """
    from astropy.io import fits
    from astropy.modeling import fitting, models
    from astropy.nddata import Cutout2D
    from astropy.wcs import WCS
    from astropy.wcs.utils import proj_plane_pixel_scales

    unit_scale = 1e6 if units == "uJy" else 1e3
    unit_label = units

    header = extract_fits_header(path)
    result: dict = {
        "Date of Observation (UTC)": header.get("DATE-OBS"),
        "Project": header.get("PROJECT"),
        "Freq (GHz)": header.get("CRVAL3"),
        f"S ({unit_label})": None,
        f"Err ({unit_label})": None,
        f"rms ({unit_label}/beam)": None,
        "SBID": header.get("SBID"),
        "Duration": header.get("DURATION"),
        "file": path,
        "ra_fit_deg": None,
        "ra_err_deg": None,
        "dec_fit_deg": None,
        "dec_err_deg": None,
        "fit_failed": False,
    }

    # ---- Load image data, squeeze degenerate axes (Stokes, frequency) ----
    with fits.open(path) as hdul:
        data_raw = hdul[0].data
        hdr = hdul[0].header

    data: np.ndarray = np.squeeze(data_raw)
    while data.ndim > 2:
        data = data[0]

    wcs2d = WCS(hdr).celestial

    # Pixel scale (arcsec/pix)
    pix_scales_deg = proj_plane_pixel_scales(wcs2d)  # [deg/pix, deg/pix]
    pixscale_arcsec = float(np.mean(np.abs(pix_scales_deg))) * 3600
    pixscale_deg = pixscale_arcsec / 3600

    # Source pixel position (column, row) = (x, y)
    x_src, y_src = (float(v) for v in wcs2d.world_to_pixel(coords))

    # ---- RMS: measure std inside an annulus around the source ----
    rms_box_pix = int(rms_radius / pixscale_arcsec * 2.0 + 1)
    try:
        rms_cutout = Cutout2D(
            data,
            position=(x_src, y_src),
            size=rms_box_pix,
            wcs=wcs2d,
            mode="partial",
            fill_value=np.nan,
        )
        rms_data = rms_cutout.data
        rx_c, ry_c = rms_cutout.input_position_cutout  # centre in cutout (x, y)
        rm_yy, rm_xx = np.mgrid[: rms_data.shape[0], : rms_data.shape[1]]
        dist_arcsec = np.sqrt((rm_xx - rx_c) ** 2 + (rm_yy - ry_c) ** 2) * pixscale_arcsec
        annulus = (dist_arcsec >= 2 * radius) & (dist_arcsec <= rms_radius) & np.isfinite(rms_data)
        rms_jybeam = float(np.std(rms_data[annulus])) if annulus.any() else float(np.nanstd(data))
    except Exception:  # noqa: BLE001
        rms_jybeam = float(np.nanstd(data))

    rms = rms_jybeam * unit_scale
    result[f"rms ({unit_label}/beam)"] = round(rms, 4)

    # ---- 2-D Gaussian fit inside circular region ----
    fit_r_pix = radius / pixscale_arcsec
    box_pix = max(int(fit_r_pix * 6 + 1), 15)  # generous box; must be odd
    if box_pix % 2 == 0:
        box_pix += 1

    try:
        cutout = Cutout2D(
            data,
            position=(x_src, y_src),
            size=box_pix,
            wcs=wcs2d,
            mode="partial",
            fill_value=np.nan,
        )
        cut = cutout.data
        cx, cy = cutout.input_position_cutout  # centre in cutout (x=col, y=row)

        good = np.isfinite(cut)
        if good.sum() < 9:
            result["fit_failed"] = True
            return result

        yy, xx = np.mgrid[: cut.shape[0], : cut.shape[1]]
        sigma_pix = fit_r_pix / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        peak_init = float(np.nanmax(cut))

        g_init = models.Gaussian2D(
            amplitude=peak_init,
            x_mean=cx,
            y_mean=cy,
            x_stddev=sigma_pix,
            y_stddev=sigma_pix,
        )

        fitter = fitting.LevMarLSQFitter(calc_uncertainties=True)
        weights = np.where(good, 1.0 / max(rms_jybeam, 1e-10), 0.0)
        g_fit = fitter(g_init, xx[good], yy[good], cut[good], weights=weights[good])

        # ierr in {1,2,3,4} means scipy considers the solution valid
        fit_failed = fitter.fit_info.get("ierr", 0) not in [1, 2, 3, 4]
        cov = fitter.fit_info.get("param_cov")  # shape (6,6): amp,x_mean,y_mean,x_std,y_std,theta

    except Exception as exc:  # noqa: BLE001
        print(f"  Fit failed: {exc}")
        result["fit_failed"] = True
        return result

    result["fit_failed"] = fit_failed
    if fit_failed:
        return result

    # ---- Extract fitted quantities ----
    peak_flux = float(g_fit.amplitude.value) * unit_scale

    # Amplitude uncertainty
    peak_err = (
        float(np.sqrt(abs(cov[0, 0]))) * unit_scale
        if cov is not None and cov.shape[0] >= 1
        else rms  # conservative fallback
    )

    combined_err = float((peak_err**2 + rms**2 + (scale * peak_flux) ** 2) ** 0.5)

    # Fitted position → sky coords
    x_fit_orig, y_fit_orig = cutout.to_original_position(
        (float(g_fit.x_mean.value), float(g_fit.y_mean.value))
    )
    sky_fit = wcs2d.pixel_to_world(x_fit_orig, y_fit_orig)
    ra_fit = float(sky_fit.ra.deg)
    dec_fit = float(sky_fit.dec.deg)

    # Position uncertainties (cov indices 1=x_mean, 2=y_mean)
    if cov is not None and cov.shape[0] >= 3:
        ra_err_deg = float(np.sqrt(abs(cov[1, 1]))) * pixscale_deg
        dec_err_deg = float(np.sqrt(abs(cov[2, 2]))) * pixscale_deg
    else:
        snr = peak_flux / rms if rms > 0 else 1.0
        pos_err_pix = fit_r_pix / (2.0 * max(snr, 1.0))
        ra_err_deg = pos_err_pix * pixscale_deg
        dec_err_deg = pos_err_pix * pixscale_deg

    result.update(
        {
            f"S ({unit_label})": round(peak_flux, 4),
            f"Err ({unit_label})": round(combined_err, 4),
            "ra_fit_deg": round(ra_fit, 6),
            "ra_err_deg": round(ra_err_deg, 6),
            "dec_fit_deg": round(dec_fit, 6),
            "dec_err_deg": round(dec_err_deg, 6),
        }
    )
    return result


# ---------------------------------------------------------------------------
# Public API — dispatches to CASA or astropy fallback
# ---------------------------------------------------------------------------


def fit_image(
    path: str,
    coords: SkyCoord,
    stokes: str = "I",
    radius: float = 5.0,
    rms_radius: float = 30.0,
    scale: float = 0.1,
    units: str = "mJy",
    workdir: str | None = None,
) -> dict:
    """
    Run forced-fit photometry on a single FITS image at the given sky position.

    Uses CASA ``imfit`` + ``imstat`` when ``casatasks`` is installed; falls back
    to a pure-Python 2-D Gaussian fit via ``astropy.modeling`` otherwise.

    Parameters
    ----------
    path : str
        Path to the FITS image.
    coords : SkyCoord
        Source coordinates.
    stokes : str
        Stokes parameter (used for RMS region in the CASA path).
    radius : float
        Fit circle radius in arcseconds.
    rms_radius : float
        Annulus outer radius in arcseconds for RMS estimation.
    scale : float
        Fractional flux-scale uncertainty added in quadrature.
    units : str
        Output flux units: ``'mJy'`` or ``'uJy'``.
    workdir : str, optional
        Directory to write CASA summary/log files (CASA path only).

    Returns
    -------
    dict
        Result row with header metadata, fitted flux, errors, and fit position.
    """
    if _CASA_AVAILABLE:
        return _fit_image_casa(path, coords, stokes, radius, rms_radius, scale, units, workdir)
    return _fit_image_astropy(path, coords, stokes, radius, rms_radius, scale, units, workdir)


# ---------------------------------------------------------------------------
# High-level entry points
# ---------------------------------------------------------------------------


def run_imfit_single(
    image_path: str,
    coords_str: str,
    output_root: str,
    name: str,
    stokes: str = "I",
    radius: float = 5.0,
    rms_radius: float = 30.0,
    scale: float = 0.1,
    units: str = "mJy",
    galactic: bool = False,
) -> dict:
    """
    Fit a single image at the given coordinates and print the result.

    Parameters
    ----------
    image_path : str
        Path to the FITS image.
    coords_str : str
        Coordinate string (e.g. ``"319.711 -57.645"``).
    output_root : str
        Root output directory.
    name : str
        Source label (used for the work subdirectory and output filename).
    stokes, radius, rms_radius, scale, units, galactic :
        Forwarded to ``fit_image()``.

    Returns
    -------
    dict
        Single-row result dict.
    """
    imfit_dir = os.path.join(output_root, name, "imfit")
    ensure_dir(imfit_dir)

    pos_eq, _ = unicoord(coords_str, galactic=galactic)
    result = fit_image(
        image_path,
        pos_eq,
        stokes=stokes,
        radius=radius,
        rms_radius=rms_radius,
        scale=scale,
        units=units,
        workdir=imfit_dir,
    )
    for k, v in result.items():
        print(f"  {k}: {v}")
    return result


def plot_imfit_lightcurve(
    csv_path: str,
    name: str,
    outdir: str,
    units: str = "mJy",
    uplim_sigma: float = 3.0,
) -> str:
    """
    Save a lightcurve plot from an imfit results CSV.

    Detections (fit_failed=False) are shown as scatter points coloured by
    frequency with error bars.  Non-detections (fit_failed=True) are shown
    as downward-pointing triangles at ``uplim_sigma`` × RMS.

    Parameters
    ----------
    csv_path : str
        Path to the CSV produced by ``run_imfit_batch()``.
    name : str
        Source label (used in title and filename).
    outdir : str
        Directory to write the plot into.
    units : str
        Flux units label (e.g. ``'mJy'``).
    uplim_sigma : float
        Upper-limit level in multiples of the local RMS.

    Returns
    -------
    str
        Path to the saved PDF, or ``''`` if there is nothing to plot.
    """
    df = pd.read_csv(csv_path)
    if df.empty:
        print("No rows in imfit CSV — nothing to plot.")
        return ""

    flux_col = f"S ({units})"
    err_col = f"Err ({units})"
    rms_col = f"rms ({units}/beam)"

    # Parse dates → MJD
    def _to_mjd(date_str: str) -> float | None:
        try:
            return float(Time(str(date_str), format="isot", scale="utc").mjd)
        except Exception:  # noqa: BLE001
            try:
                return float(Time(str(date_str)).mjd)
            except Exception:  # noqa: BLE001
                return None

    df["mjd"] = df["Date of Observation (UTC)"].apply(_to_mjd)
    df = df.dropna(subset=["mjd"]).sort_values("mjd").reset_index(drop=True)

    if df.empty:
        print("Could not parse any observation dates — skipping plot.")
        return ""

    # Frequency in GHz for colouring
    freq_col = "Freq (GHz)"
    has_freq = freq_col in df.columns and df[freq_col].notna().any()
    if has_freq:
        norm = plt.Normalize(df[freq_col].min(), df[freq_col].max())

    detected = df[df[flux_col].notna()].copy()
    non_detected = df[df[flux_col].isna()].copy()

    fig, ax = plt.subplots()

    # ---- Detections ----
    if not detected.empty:
        det = detected
        sc = ax.scatter(
            det["mjd"],
            det[flux_col],
            c=det[freq_col] if has_freq else None,
            cmap="viridis_r" if has_freq else None,
            norm=norm if has_freq else None,
            s=50,
            zorder=3,
        )
        if err_col in det.columns:
            ax.errorbar(
                det["mjd"],
                det[flux_col],
                yerr=det[err_col].fillna(0),
                fmt="none",
                ecolor="gray",
                capsize=2,
                alpha=0.5,
                zorder=2,
            )
        if has_freq:
            cbar = fig.colorbar(sc, ax=ax)
            cbar.set_label("Frequency [GHz]")

    # ---- Upper limits (rows with no flux measurement: 3σ × RMS) ----
    if not non_detected.empty and rms_col in non_detected.columns:
        ul = non_detected.dropna(subset=[rms_col])
        if not ul.empty:
            ax.scatter(
                ul["mjd"],
                ul[rms_col] * uplim_sigma,
                marker="v",
                s=60,
                c="gray",
                alpha=0.7,
                zorder=3,
                label=f"{uplim_sigma:.0f}σ upper limits",
            )
            ax.legend(loc="upper right", fontsize="small")

    ax.set_xlabel("Date [MJD]")
    ax.set_ylabel(f"Flux Density [{units}/beam]")
    ax.set_title(name)
    fig.tight_layout()

    timestamp = datetime.now().strftime("%Y%m%d_%H%M")
    save_path = os.path.join(outdir, f"{name}_imfit_lc_{timestamp}.pdf")
    fig.savefig(save_path, dpi=300)
    plt.close(fig)
    print(f"Lightcurve plot saved: {save_path}")
    return save_path


def run_imfit_batch(
    image_paths: list[str],
    coords_str: str,
    output_root: str,
    name: str,
    stokes: str = "I",
    radius: float = 5.0,
    rms_radius: float = 30.0,
    scale: float = 0.1,
    units: str = "mJy",
    galactic: bool = False,
    plot: bool = True,
    pmra: float = 0.0,
    pmdec: float = 0.0,
    epoch: str = "J2000",
) -> str:
    """
    Fit multiple images at the same sky position and write results to CSV.

    Parameters
    ----------
    image_paths : list[str]
        Paths to FITS images (one row per image in the output CSV).
    coords_str : str
        Coordinate string.
    output_root : str
        Root output directory.
    name : str
        Source label.
    stokes, radius, rms_radius, scale, units, galactic :
        Forwarded to ``fit_image()``.
    plot : bool
        If ``True`` (default), save a lightcurve plot after writing the CSV.
    pmra, pmdec : float
        Proper motion in mas/yr.  Set both to 0 to disable.
    epoch : str
        Reference epoch for the catalogue position (e.g. ``"J2000"``).

    Returns
    -------
    str
        Path to the output CSV.
    """
    backend = "CASA" if _CASA_AVAILABLE else "astropy"
    print(f"Using {backend} backend for imfit.")

    imfit_dir = os.path.join(output_root, name, "imfit")
    ensure_dir(imfit_dir)

    pos_eq, _ = unicoord(coords_str, galactic=galactic)
    apply_pm = pmra != 0.0 or pmdec != 0.0
    rows = []
    for i, path in enumerate(image_paths, 1):
        print(f"[{i}/{len(image_paths)}] Fitting {Path(path).name}…")
        # Compute PM-corrected position for this image's epoch
        coords = pos_eq
        if apply_pm:
            try:
                hdr = extract_fits_header(path)
                mjd_obs = hdr.get("MJD-OBS")
                if mjd_obs is None:
                    date_obs = hdr.get("DATE-OBS")
                    if date_obs:
                        mjd_obs = float(Time(str(date_obs)).mjd)
                if mjd_obs is not None:
                    coords = apply_proper_motion(pos_eq, pmra, pmdec, epoch, float(mjd_obs))
            except Exception:  # noqa: BLE001
                coords = pos_eq
        rows.append(
            fit_image(
                path,
                coords,
                stokes=stokes,
                radius=radius,
                rms_radius=rms_radius,
                scale=scale,
                units=units,
                workdir=imfit_dir,
            )
        )

    if not rows:
        print("No images processed.")
        return ""

    csv_path = os.path.join(imfit_dir, f"{name}_imfit_results.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nSaved {len(rows)} rows to {csv_path}")

    if plot:
        plot_imfit_lightcurve(csv_path, name, imfit_dir, units=units)

    return csv_path
