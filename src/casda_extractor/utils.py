"""Shared utility functions for coordinate handling and FITS header extraction."""

import math
import re
import subprocess

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time


def unicoord(coords: str, galactic: bool = False, display: bool = False) -> tuple:
    """
    Universal coordinate converter for equatorial and galactic coordinates.

    Parameters
    ----------
    coords : str
        Source coordinates in any standard format (e.g. "10.684 41.268" or "00:42:44 +41:16:09")
    galactic : bool
        If True, assumes galactic coordinates (l, b).
    display : bool
        If True, print all coordinate representations.

    Returns
    -------
    tuple
        (pos_eq, pos_gal) as SkyCoord objects.
    """
    ra, dec = coords.split(" ")
    raunit = "hourangle" if (":" in ra or "h" in ra) else "deg"
    if not galactic:
        pos_eq = SkyCoord(ra=ra, dec=dec, unit=(raunit, "deg"))
        pos_gal = pos_eq.galactic
    else:
        pos_gal = SkyCoord(l=ra, b=dec, unit=("deg", "deg"), frame="galactic")
        pos_eq = pos_gal.transform_to("icrs")

    if display:
        print(29 * "=", "\nEquatorial Coordinates:")
        print(pos_eq.to_string(style="decimal", precision=6))
        print(pos_eq.to_string(style="hmsdms", precision=3))
        print(29 * "=", "\nGalactic Coordinates:")
        print(pos_gal.to_string(style="decimal", precision=6))
        print(29 * "=")
    return pos_eq, pos_gal


def extract_fits_header(fits_path: str) -> dict:
    """
    Extract key header values from a FITS file.

    Uses ``fitshdr`` if available, falls back to astropy's ``fitsheader``.

    Parameters
    ----------
    fits_path : str
        Path to the FITS file.

    Returns
    -------
    dict
        Keys: DATE-OBS, PROJECT, CRVAL3 (GHz), SBID, DURATION (minutes).
    """
    keys = ["DATE-OBS", "PROJECT", "CRVAL3", "SBID", "DURATION"]

    try:
        result = subprocess.run(
            ["fitshdr", fits_path],
            capture_output=True,
            text=True,
            check=True,
        )
        header_text = result.stdout
    except FileNotFoundError:
        result = subprocess.run(
            ["fitsheader", fits_path],
            capture_output=True,
            text=True,
            check=True,
        )
        header_text = result.stdout
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Failed to read FITS header: {e.stderr}") from e

    values: dict = {}
    for key in keys:
        pattern = rf"^{key}\s*=\s*([^/\n]+)"
        match = re.search(pattern, header_text, re.MULTILINE)
        if match:
            raw = match.group(1).strip().strip("'").strip()
            try:
                values[key] = float(raw) if ("." in raw or "e" in raw.lower()) else int(raw)
            except ValueError:
                values[key] = raw
        else:
            values[key] = None

    if values.get("CRVAL3") is not None:
        values["CRVAL3"] = values["CRVAL3"] / 1e9  # Hz → GHz

    if values.get("DURATION") is not None:
        values["DURATION"] = values["DURATION"] / 60  # s → min

    return values


def apply_proper_motion(
    coord: SkyCoord,
    pmra: float,
    pmdec: float,
    epoch: str,
    obs_mjd: float,
) -> SkyCoord:
    """
    Return ``coord`` propagated to ``obs_mjd`` via proper motion.

    Falls back to the original ``coord`` on any error (including non-finite
    ``obs_mjd``).

    Parameters
    ----------
    coord : SkyCoord
        Base ICRS position.
    pmra : float
        Proper motion in RA × cos δ (mas/yr).
    pmdec : float
        Proper motion in Dec (mas/yr).
    epoch : str
        Reference epoch string (e.g. ``"J2000"`` or ``"2016.0"``).
    obs_mjd : float
        Observation epoch in MJD.

    Returns
    -------
    SkyCoord
        PM-corrected position, or ``coord`` if the correction fails.
    """
    try:
        if not math.isfinite(obs_mjd):
            return coord
        return SkyCoord(
            coord.ra.deg * u.deg,
            coord.dec.deg * u.deg,
            pm_ra_cosdec=pmra * u.mas / u.yr,
            pm_dec=pmdec * u.mas / u.yr,
            frame="icrs",
            obstime=Time(epoch),
            distance=1.0 * u.pc,
        ).apply_space_motion(Time(obs_mjd, format="mjd"))
    except Exception:  # noqa: BLE001
        return coord


def parse_coords_input(coords_str: str, galactic: bool = False) -> SkyCoord:
    """
    Parse a coordinate string into an ICRS SkyCoord.

    Parameters
    ----------
    coords_str : str
        e.g. "10.684 41.268" or "00:42:44 +41:16:09"
    galactic : bool
        Treat input as galactic (l, b) coordinates.

    Returns
    -------
    SkyCoord
        ICRS position.
    """
    pos_eq, _ = unicoord(coords_str, galactic=galactic)
    return pos_eq
