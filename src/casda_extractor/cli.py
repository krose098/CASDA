"""
casda — CLI toolkit for CASDA archive queries and forced photometry.

Subcommands
-----------
lightcurves   Query and download detection lightcurves from CASDA catalogues.
images        Download full-field images and generate source cutouts.
imfit         Forced-fit photometry on FITS images using CASA imfit/imstat.
cleanup       Remove large intermediate files while keeping results.

Single-source usage (--coords / -c):
    casda lightcurves -c "10.684 41.268" -n M31 -o ./results

Multi-source usage (--csv):
    casda images --csv sources.csv -o ./results

sources.csv must contain at minimum ``ra`` and ``dec`` columns;
an optional ``name`` column is used as the source label.
"""

from __future__ import annotations

import os
import sys

import click
import pandas as pd

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


def _load_sources_csv(csv_path: str) -> pd.DataFrame:
    """Read a sources CSV and validate required columns."""
    df = pd.read_csv(csv_path)
    for col in ("ra", "dec"):
        if col not in df.columns:
            click.echo(f"ERROR: CSV is missing required column '{col}'.", err=True)
            sys.exit(1)
    if "name" not in df.columns:
        # Auto-generate names from coordinates
        df["name"] = df.apply(lambda r: f"J{r['ra']:.3f}{'+' if r['dec'] >= 0 else ''}{r['dec']:.3f}", axis=1)
        click.echo("No 'name' column found — auto-generating source names from coordinates.")
    return df


def _single_source_df(coords: str, name: str, galactic: bool) -> tuple[float, float, str]:
    """Parse --coords into (ra_deg, dec_deg, name)."""
    from .utils import unicoord

    pos_eq, _ = unicoord(coords, galactic=galactic)
    ra = float(pos_eq.ra.deg)
    dec = float(pos_eq.dec.deg)
    if not name:
        name = f"J{ra:.3f}{'+' if dec >= 0 else ''}{dec:.3f}"
    return ra, dec, name


# ---------------------------------------------------------------------------
# Root group
# ---------------------------------------------------------------------------


@click.group()
@click.version_option()
def cli() -> None:
    """CASDA CLI toolkit for radio transient and variable star studies."""


# ---------------------------------------------------------------------------
# lightcurves
# ---------------------------------------------------------------------------


@cli.command("lightcurves")
@click.option("-c", "--coords", default=None, help="Source coordinates: 'RA DEC' in degrees (or HH:MM:SS).")
@click.option("-n", "--name", default="", help="Source name / label.")
@click.option(
    "--csv",
    "csv_path",
    default=None,
    type=click.Path(exists=True),
    help="CSV file with columns: ra, dec[, name, ...].",
)
@click.option(
    "-u", "--username", default=None, envvar="CASDA_USER", help="OPAL username (or set CASDA_USER env var)."
)
@click.option("-o", "--output", default="./casda_output", show_default=True, help="Root output directory.")
@click.option("--search-radius", default=10.0, show_default=True, help="TAP cone search radius in degrees.")
@click.option("--sep", default=5.0, show_default=True, help="Cross-match radius in arcseconds.")
@click.option(
    "--scale", default=0.06, show_default=True, help="Fractional flux-scale uncertainty for quadrature error."
)
@click.option("--max-workers", default=8, show_default=True, help="Parallel download threads.")
@click.option("--no-plot", is_flag=True, default=False, help="Skip lightcurve plot.")
@click.option("-g", "--galactic", is_flag=True, default=False, help="Treat --coords as galactic (l, b).")
@click.option("--pmra", default=0.0, show_default=True, help="Proper motion in RA × cos δ (mas/yr).")
@click.option("--pmdec", default=0.0, show_default=True, help="Proper motion in Dec (mas/yr).")
@click.option("--epoch", default="J2000", show_default=True, help="Reference epoch for the source position.")
def cmd_lightcurves(
    coords,
    name,
    csv_path,
    username,
    output,
    search_radius,
    sep,
    scale,
    max_workers,
    no_plot,
    galactic,
    pmra,
    pmdec,
    epoch,
):
    """Query CASDA for detection lightcurves (catalogue cross-matching)."""
    from .auth import get_casda
    from .lightcurves import run_lightcurves

    casda_obj = get_casda(username)

    if csv_path:
        df = _load_sources_csv(csv_path)
        for _, row in df.iterrows():
            run_lightcurves(
                name=str(row["name"]),
                ra=float(row["ra"]),
                dec=float(row["dec"]),
                output_root=output,
                casda_obj=casda_obj,
                search_radius=search_radius,
                sep=sep,
                scale=scale,
                max_workers=max_workers,
                plot=not no_plot,
                pmra=float(row.get("pmra", 0.0) or 0.0),
                pmdec=float(row.get("pmdec", 0.0) or 0.0),
                epoch=str(row.get("epoch", "J2000") or "J2000"),
            )
    elif coords:
        ra, dec, src_name = _single_source_df(coords, name, galactic)
        run_lightcurves(
            name=src_name,
            ra=ra,
            dec=dec,
            output_root=output,
            casda_obj=casda_obj,
            search_radius=search_radius,
            sep=sep,
            scale=scale,
            max_workers=max_workers,
            plot=not no_plot,
            pmra=pmra,
            pmdec=pmdec,
            epoch=epoch,
        )
    else:
        click.echo("Provide either --coords or --csv.", err=True)
        sys.exit(1)


# ---------------------------------------------------------------------------
# images
# ---------------------------------------------------------------------------


@cli.command("images")
@click.option("-c", "--coords", default=None, help="Source coordinates: 'RA DEC' in degrees.")
@click.option("-n", "--name", default="", help="Source name / label.")
@click.option(
    "--csv",
    "csv_path",
    default=None,
    type=click.Path(exists=True),
    help="CSV file with columns: ra, dec[, name, pmra, pmdec, epoch].",
)
@click.option(
    "-u", "--username", default=None, envvar="CASDA_USER", help="OPAL username (or set CASDA_USER env var)."
)
@click.option("-o", "--output", default="./casda_output", show_default=True, help="Root output directory.")
@click.option("--search-radius", default=10.0, show_default=True, help="Image search cone in degrees.")
@click.option("--cutout-size", default=4.0, show_default=True, help="Cutout size in arcminutes.")
@click.option("--stokes", default="I", show_default=True, help="Stokes parameter.")
@click.option("--max-workers", default=8, show_default=True, help="Parallel download threads.")
@click.option("--plot-only", is_flag=True, default=False, help="Skip download and re-plot existing cutouts.")
@click.option("-g", "--galactic", is_flag=True, default=False, help="Treat --coords as galactic (l, b).")
@click.option("--pmra", default=0.0, show_default=True, help="Proper motion in RA × cos δ (mas/yr).")
@click.option("--pmdec", default=0.0, show_default=True, help="Proper motion in Dec (mas/yr).")
@click.option("--epoch", default="J2000", show_default=True, help="Reference epoch for the source position.")
def cmd_images(
    coords,
    name,
    csv_path,
    username,
    output,
    search_radius,
    cutout_size,
    stokes,
    max_workers,
    plot_only,
    galactic,
    pmra,
    pmdec,
    epoch,
):
    """Download ASKAP images and generate source cutouts."""
    from .auth import get_casda
    from .images import run_images

    casda_obj = None if plot_only else get_casda(username)

    if csv_path:
        df = _load_sources_csv(csv_path)
        for _, row in df.iterrows():
            run_images(
                name=str(row["name"]),
                ra=float(row["ra"]),
                dec=float(row["dec"]),
                output_root=output,
                casda_obj=casda_obj,
                pmra=float(row.get("pmra", 0.0) or 0.0),
                pmdec=float(row.get("pmdec", 0.0) or 0.0),
                epoch=str(row.get("epoch", "J2000") or "J2000"),
                stokes=stokes,
                search_radius=search_radius,
                cutout_size=cutout_size,
                max_workers=max_workers,
                plot_only=plot_only,
            )
    elif coords:
        ra, dec, src_name = _single_source_df(coords, name, galactic)
        run_images(
            name=src_name,
            ra=ra,
            dec=dec,
            output_root=output,
            casda_obj=casda_obj,
            pmra=pmra,
            pmdec=pmdec,
            epoch=epoch,
            stokes=stokes,
            search_radius=search_radius,
            cutout_size=cutout_size,
            max_workers=max_workers,
            plot_only=plot_only,
        )
    else:
        click.echo("Provide either --coords or --csv.", err=True)
        sys.exit(1)


# ---------------------------------------------------------------------------
# imfit
# ---------------------------------------------------------------------------


@cli.group("imfit", invoke_without_command=True)
@click.pass_context
@click.option(
    "-c", "--coords", default=None, help="Source coordinates: 'RA DEC' in degrees. Required unless --csv."
)
@click.option("-n", "--name", default="", help="Source name / label.")
@click.option(
    "--csv",
    "csv_path",
    default=None,
    type=click.Path(exists=True),
    help="CSV with columns: ra, dec[, name]. Runs batch mode per source.",
)
@click.option("-o", "--output", default="./casda_output", show_default=True, help="Root output directory.")
@click.option(
    "--images-dir", default=None, help="Directory containing FITS images (overrides output/<name>/images/)."
)
@click.option("--stokes", default="I", show_default=True, help="Stokes parameter.")
@click.option("--radius", default=5.0, show_default=True, help="Fit circle radius in arcseconds.")
@click.option(
    "--rms-radius", default=30.0, show_default=True, help="Annulus outer radius for RMS in arcseconds."
)
@click.option("--scale", default=0.1, show_default=True, help="Fractional flux-scale uncertainty.")
@click.option(
    "-u",
    "--units",
    type=click.Choice(["mJy", "uJy"]),
    default="mJy",
    show_default=True,
    help="Output flux units.",
)
@click.option("--no-plot", is_flag=True, default=False, help="Skip lightcurve plot.")
@click.option("-g", "--galactic", is_flag=True, default=False, help="Treat --coords as galactic (l, b).")
@click.option("--pmra", default=0.0, show_default=True, help="Proper motion in RA × cos δ (mas/yr).")
@click.option("--pmdec", default=0.0, show_default=True, help="Proper motion in Dec (mas/yr).")
@click.option("--epoch", default="J2000", show_default=True, help="Reference epoch for the source position.")
def cmd_imfit(
    ctx,
    coords,
    name,
    csv_path,
    output,
    images_dir,
    stokes,
    radius,
    rms_radius,
    scale,
    units,
    no_plot,
    galactic,
    pmra,
    pmdec,
    epoch,
):
    """
    Forced-fit photometry on downloaded FITS images.

    Images are located under OUTPUT/NAME/images/ (or --images-dir).

    \b
    Backends (selected automatically):
      CASA      casatasks.imfit + imstat — highest fidelity; used when
                casatasks is importable (Linux, pip install casda[casa]).
      astropy   astropy.modeling Gaussian2D — pure Python, works on all
                platforms including macOS; installed by default.

    Both backends produce the same output columns and CSV schema.
    The active backend is reported at the start of each run.
    """
    if ctx.invoked_subcommand is not None:
        return
    import glob

    from .imfit import run_imfit_batch

    def _get_images(src_name: str) -> list[str]:
        img_root = images_dir or os.path.join(output, src_name, "images")
        patterns = ["*.fits", "*.FITS"]
        files: list[str] = []
        for pat in patterns:
            files += glob.glob(os.path.join(img_root, pat))
        files.sort()
        if not files:
            click.echo(f"No FITS images found under {img_root}", err=True)
        return files

    if csv_path:
        df = _load_sources_csv(csv_path)
        for _, row in df.iterrows():
            src_name = str(row["name"])
            coords_str = f"{row['ra']} {row['dec']}"
            files = _get_images(src_name)
            if files:
                run_imfit_batch(
                    image_paths=files,
                    coords_str=coords_str,
                    output_root=output,
                    name=src_name,
                    stokes=stokes,
                    radius=radius,
                    rms_radius=rms_radius,
                    scale=scale,
                    units=units,
                    plot=not no_plot,
                    pmra=float(row.get("pmra", 0.0) or 0.0),
                    pmdec=float(row.get("pmdec", 0.0) or 0.0),
                    epoch=str(row.get("epoch", "J2000") or "J2000"),
                )
    elif coords:
        from .utils import unicoord

        pos_eq, _ = unicoord(coords, galactic=galactic)
        ra_deg = float(pos_eq.ra.deg)
        dec_deg = float(pos_eq.dec.deg)
        src_name = name or f"J{ra_deg:.3f}{'+' if dec_deg >= 0 else ''}{dec_deg:.3f}"
        coords_str = f"{ra_deg} {dec_deg}"
        files = _get_images(src_name)
        if files:
            run_imfit_batch(
                image_paths=files,
                coords_str=coords_str,
                output_root=output,
                name=src_name,
                stokes=stokes,
                radius=radius,
                rms_radius=rms_radius,
                scale=scale,
                units=units,
                plot=not no_plot,
                pmra=pmra,
                pmdec=pmdec,
                epoch=epoch,
            )
    else:
        click.echo("Provide either --coords or --csv.", err=True)
        sys.exit(1)


@cmd_imfit.command("plot")
@click.option(
    "--csv",
    "csv_path",
    required=True,
    type=click.Path(exists=True),
    help="Path to an imfit results CSV.",
)
@click.option("-n", "--name", default="", help="Source name (used in title and filename).")
@click.option("-o", "--output", default=None, help="Output directory (defaults to CSV directory).")
@click.option(
    "-u",
    "--units",
    type=click.Choice(["mJy", "uJy"]),
    default="mJy",
    show_default=True,
    help="Flux units label (must match the CSV).",
)
@click.option(
    "--uplim-sigma",
    default=3.0,
    show_default=True,
    help="Upper-limit level in multiples of the local RMS.",
)
def cmd_imfit_plot(csv_path, name, output, units, uplim_sigma):
    """Plot a lightcurve from an existing imfit results CSV.

    Rows with a flux measurement are shown as scatter points (coloured by
    frequency when available).  Rows with no flux measurement are shown as
    downward-pointing triangles at UPLIM_SIGMA × RMS.
    """
    from .imfit import plot_imfit_lightcurve

    outdir = output or os.path.dirname(os.path.abspath(csv_path))
    src_name = name or os.path.splitext(os.path.basename(csv_path))[0]
    plot_imfit_lightcurve(csv_path, src_name, outdir, units=units, uplim_sigma=uplim_sigma)


# ---------------------------------------------------------------------------
# credentials
# ---------------------------------------------------------------------------


@cli.group("credentials")
def cmd_credentials() -> None:
    """Manage stored OPAL credentials in the system keyring."""


@cmd_credentials.command("save")
@click.option("-u", "--username", prompt="OPAL username (email)", help="OPAL username.")
@click.option(
    "-p", "--password", prompt=True, hide_input=True, confirmation_prompt=True, help="OPAL password."
)
def credentials_save(username: str, password: str) -> None:
    """Save OPAL credentials to the system keyring."""
    from .auth import save_credentials

    save_credentials(username, password)


@cmd_credentials.command("delete")
def credentials_delete() -> None:
    """Remove stored OPAL credentials from the system keyring."""
    from .auth import delete_credentials

    delete_credentials()


@cmd_credentials.command("show")
def credentials_show() -> None:
    """Show the username currently stored in the system keyring (password is hidden)."""
    from .auth import load_credentials

    username, password = load_credentials()
    if username is None:
        click.echo("No credentials stored.")
    else:
        click.echo(f"Stored username: {username}")
        click.echo("Password: <hidden>")


# ---------------------------------------------------------------------------
# cleanup
# ---------------------------------------------------------------------------


@cli.command("cleanup")
@click.option("-n", "--name", default=None, help="Source name to clean up.  Omit to clean all sources.")
@click.option("-o", "--output", default="./casda_output", show_default=True, help="Root output directory.")
@click.option(
    "--dry-run", is_flag=True, default=False, help="Print what would be deleted without removing anything."
)
@click.option(
    "--keep-images", is_flag=True, default=False, help="Keep full-field FITS images (only remove catalogues)."
)
@click.option(
    "--keep-catalogues", is_flag=True, default=False, help="Keep catalogue XML files (only remove images)."
)
def cmd_cleanup(name, output, dry_run, keep_images, keep_catalogues):
    """
    Remove large intermediate files (images/ and/or catalogues/) for one or all sources.

    Lightcurve CSVs, cutouts, imfit results, and plots are always preserved.
    """
    from .output import cleanup_all, cleanup_source

    remove_images = not keep_images
    remove_catalogues = not keep_catalogues

    if name:
        cleanup_source(
            output,
            name,
            dry_run=dry_run,
            remove_images=remove_images,
            remove_catalogues=remove_catalogues,
        )
    else:
        cleanup_all(
            output,
            dry_run=dry_run,
            remove_images=remove_images,
            remove_catalogues=remove_catalogues,
        )


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> None:
    """Package entry point."""
    cli()


if __name__ == "__main__":
    main()
