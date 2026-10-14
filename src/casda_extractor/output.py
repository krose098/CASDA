"""Output directory management and cleanup utilities."""

from __future__ import annotations

import contextlib
import os
import shutil


def ensure_dir(path: str) -> str:
    """Create ``path`` if it does not exist. Returns the path."""
    os.makedirs(path, exist_ok=True)
    return path


def source_dir(output_root: str, name: str) -> str:
    """Return (and create) the per-source output directory."""
    return ensure_dir(os.path.join(output_root, name))


def _human_size(nbytes: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if nbytes < 1024:
            return f"{nbytes:.1f} {unit}"
        nbytes //= 1024
    return f"{nbytes:.1f} TB"


def _dir_size(path: str) -> int:
    total = 0
    for root, _, files in os.walk(path):
        for f in files:
            with contextlib.suppress(OSError):
                total += os.path.getsize(os.path.join(root, f))
    return total


# ---------------------------------------------------------------------------
# Cleanup
# ---------------------------------------------------------------------------


SECONDARY_SUBDIRS = ("images", "catalogues")
"""Subdirectories that hold large downloaded data (safe to delete after processing)."""


def cleanup_source(
    output_root: str,
    name: str,
    dry_run: bool = False,
    remove_images: bool = True,
    remove_catalogues: bool = True,
) -> int:
    """
    Delete secondary/intermediate data for one source, keeping results and cutouts.

    Parameters
    ----------
    output_root : str
        Root output directory.
    name : str
        Source label.
    dry_run : bool
        If True, print what would be deleted without removing anything.
    remove_images : bool
        Delete the ``images/`` subdirectory (full-field FITS files).
    remove_catalogues : bool
        Delete the ``catalogues/`` subdirectory (downloaded XML files).

    Returns
    -------
    int
        Bytes freed (0 in dry-run mode).
    """
    src_dir = os.path.join(output_root, name)
    if not os.path.isdir(src_dir):
        print(f"No output directory found for '{name}' at {src_dir}")
        return 0

    targets = []
    if remove_images:
        targets.append("images")
    if remove_catalogues:
        targets.append("catalogues")

    freed = 0
    for subdir in targets:
        target = os.path.join(src_dir, subdir)
        if os.path.isdir(target):
            size = _dir_size(target)
            if dry_run:
                print(f"  [dry-run] Would remove: {target}  ({_human_size(size)})")
            else:
                shutil.rmtree(target)
                freed += size
                print(f"  Removed: {target}  ({_human_size(size)} freed)")
        else:
            print(f"  Not found (already clean?): {target}")

    if not dry_run and freed:
        print(f"Freed {_human_size(freed)} for source '{name}'.")
    return freed


def cleanup_all(
    output_root: str,
    dry_run: bool = False,
    remove_images: bool = True,
    remove_catalogues: bool = True,
) -> int:
    """
    Run ``cleanup_source`` for every subdirectory found under ``output_root``.

    Parameters
    ----------
    output_root : str
        Root output directory containing one subdirectory per source.
    dry_run : bool
        Print only, do not delete.
    remove_images, remove_catalogues : bool
        Which secondary subdirectory types to remove.

    Returns
    -------
    int
        Total bytes freed.
    """
    if not os.path.isdir(output_root):
        print(f"Output root not found: {output_root}")
        return 0

    names = [d for d in os.listdir(output_root) if os.path.isdir(os.path.join(output_root, d))]
    if not names:
        print("No source directories found.")
        return 0

    total = 0
    for name in sorted(names):
        print(f"\n[{name}]")
        total += cleanup_source(
            output_root,
            name,
            dry_run=dry_run,
            remove_images=remove_images,
            remove_catalogues=remove_catalogues,
        )

    if total:
        print(f"\nTotal freed: {_human_size(total)}")
    return total
