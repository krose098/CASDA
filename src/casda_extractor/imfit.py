"""Forced-fit photometry using CASA's imfit / imstat tasks."""

from __future__ import annotations

import csv
import os
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord

from .output import ensure_dir
from .utils import extract_fits_header, unicoord

try:
    from casatasks import imfit, imstat  # type: ignore[import]

    _CASA_AVAILABLE = True
except ImportError:
    _CASA_AVAILABLE = False


def _check_casa() -> None:
    if not _CASA_AVAILABLE:
        raise RuntimeError(
            "casatasks is not installed or could not be imported. "
            "Install CASA (https://casa.nrao.edu/) and make sure casatasks is on your PYTHONPATH."
        )


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


# ---------------------------------------------------------------------------
# Core extraction
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
    Run CASA imfit + imstat on a single FITS image at the given sky position.

    Parameters
    ----------
    path : str
        Path to the FITS image.
    coords : SkyCoord
        Source coordinates.
    stokes : str
        Stokes parameter (for imstat RMS measurement).
    radius : float
        Fit circle radius in arcseconds.
    rms_radius : float
        Annulus outer radius in arcseconds for RMS estimation.
    scale : float
        Fractional flux-scale uncertainty added in quadrature.
    units : str
        Output flux units: ``'mJy'`` or ``'uJy'``.
    workdir : str, optional
        Directory to write CASA summary/log files.  Defaults to the current directory.

    Returns
    -------
    dict
        Result row with header metadata, fitted flux, errors, and fit position.
    """
    _check_casa()
    unit_scale = 1e6 if units == "uJy" else 1e3
    unit_label = units  # 'mJy' or 'uJy'

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

    # Derive a short source label for CASA output files
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

    Returns
    -------
    str
        Path to the output CSV.
    """
    imfit_dir = os.path.join(output_root, name, "imfit")
    ensure_dir(imfit_dir)

    pos_eq, _ = unicoord(coords_str, galactic=galactic)
    rows = []
    for i, path in enumerate(image_paths, 1):
        print(f"[{i}/{len(image_paths)}] Fitting {Path(path).name}…")
        rows.append(
            fit_image(
                path,
                pos_eq,
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
    return csv_path
