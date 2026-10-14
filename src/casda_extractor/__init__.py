"""
casda_extractor — CASDA querying and forced-photometry toolkit.

Public API
----------
auth          CASDA session management.
lightcurves   Detection lightcurve extraction.
images        Image download and cutout generation.
imfit         Forced-fit photometry via CASA.
output        Output directory management and cleanup.
utils         Coordinate parsing and FITS header utilities.
"""

from __future__ import annotations

import importlib

__all__ = [
    "get_casda",
    "reset_casda",
    "ensure_dir",
    "cleanup_source",
    "cleanup_all",
    "unicoord",
    "parse_coords_input",
    "extract_fits_header",
]

_LAZY: dict[str, tuple[str, str]] = {
    "get_casda": (".auth", "get_casda"),
    "reset_casda": (".auth", "reset_casda"),
    "ensure_dir": (".output", "ensure_dir"),
    "cleanup_source": (".output", "cleanup_source"),
    "cleanup_all": (".output", "cleanup_all"),
    "unicoord": (".utils", "unicoord"),
    "parse_coords_input": (".utils", "parse_coords_input"),
    "extract_fits_header": (".utils", "extract_fits_header"),
}


def __getattr__(name: str) -> object:
    if name in _LAZY:
        mod_name, attr = _LAZY[name]
        mod = importlib.import_module(mod_name, package=__name__)
        return getattr(mod, attr)
    raise AttributeError(f"module 'casda_extractor' has no attribute {name!r}")
