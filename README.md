# CASDA
## a repository for CASDA-related tools

I originally developed this notebook (based on code obtained from Dr. Laura Driessen, who got it from Dr. Minh Huynh) to download publicly available CASDA data around a given object's coordinates (or a .csv of objects/coordinates). [View Notebook](https://github.com/krose098/CASDA/blob/main/CASDA_Lightcurve_Public.ipynb). Below see the installation and usage guide for the command line interface (CLI) version.

All you need to do is have an OPAL login, add your source details, (pray that CASDA is working), and it will:
1) Download the component catalogues into a newly created, named folder (you can play around with the code and download other data products too)
2) Match components to within a given sky separation (default 5")
3) Produce a lightcurve from all public ASKAP data
4) Save the lightcurve data points to a .csv.

A couple of small things to keep in mind with this data:
* The matches are to Selavy component source catalogues (within 5") of the object coordinates and Selavy has a 5sigma detection threshold.
* There was a recent bug fix due to a change in the CASDA file naming convention. It is possible that some files are missed by my search if they are unconventionally named.
* In some cases CASDA may store multiple versions of a single observation -- especially from pilot programs -- if the data was reprocessed. This may produce, for example, two detections on the same day with slightly different fluxes.
* The columns in each .csv should be self-explanatory (feel free to reach out if not) with the exception of flux_err_quad which is a more conservation flux uncertainty I defined for the lightcurves as: $$((flux_{err}^2) + (rms^2) + (0.06*flux)^2)^{0.5}$$

## CLI installation

```bash
pip install -e .          # core install (all platforms)
pip install -e '.[dev]'   # + dev tools
pip install -e '.[casa]'  # + CASA casatasks (Linux only; see below)
```

## CLI usage

All commands accept either a single source via `--coords` / `-c` or a batch of sources via `--csv`.
A sources CSV must contain at minimum `ra` and `dec` columns; an optional `name` column is used as the source label.

### casda credentials

Store your OPAL login in the system keyring so you don't have to type it on every run.

```bash
casda credentials save                        # prompts for username and password
casda credentials show                        # print stored username (password hidden)
casda credentials delete                      # remove stored credentials
```

You can also set the `CASDA_USER` environment variable to supply your username without the keyring.

### casda lightcurves

Query CASDA for Selavy component catalogues near a source, cross-match to within a given sky separation, and build a detection lightcurve. Results are saved as a CSV and (by default) a scatter plot coloured by frequency.

```bash
# single source
casda lightcurves -c "10.684 41.268" -n M31 -o ./results

# galactic coordinates
casda lightcurves -c "121.17 -21.57" -n CasA -o ./results --galactic

# batch from CSV (columns: ra, dec, name)
casda lightcurves --csv sources.csv -o ./results

# skip the plot
casda lightcurves -c "10.684 41.268" -n M31 -o ./results --no-plot

# widen the search cone or tighten the cross-match
casda lightcurves -c "10.684 41.268" -n M31 -o ./results --search-radius 15 --sep 3
```

Key options:

| Option | Default | Description |
|--------|---------|-------------|
| `--search-radius` | 10.0° | TAP cone search radius |
| `--sep` | 5.0″ | Cross-match radius |
| `--scale` | 0.06 | Fractional flux-scale uncertainty |
| `--max-workers` | 8 | Parallel download threads |
| `--no-plot` | off | Skip the lightcurve PDF |

Outputs land in `OUTPUT/NAME/lightcurves/`.

### casda images

Download full-field ASKAP FITS images around a source and generate postage-stamp cutouts. Proper-motion correction is applied when `pmra`/`pmdec`/`epoch` are supplied in the CSV.

```bash
# single source
casda images -c "10.684 41.268" -n M31 -o ./results

# batch from CSV (optional columns: pmra, pmdec, epoch)
casda images --csv sources.csv -o ./results

# re-generate cutout plots from already-downloaded images (no re-download)
casda images -c "10.684 41.268" -n M31 -o ./results --plot-only

# 2′ cutout, Stokes V
casda images -c "10.684 41.268" -n M31 -o ./results --cutout-size 2 --stokes V
```

Key options:

| Option | Default | Description |
|--------|---------|-------------|
| `--search-radius` | 10.0° | Image search cone |
| `--cutout-size` | 4.0′ | Postage-stamp side length |
| `--stokes` | I | Stokes parameter |
| `--max-workers` | 8 | Parallel download threads |
| `--plot-only` | off | Re-plot from existing FITS files |

Full-field images land in `OUTPUT/NAME/images/`; cutouts in `OUTPUT/NAME/cutouts/`.

### casda imfit

Fit a 2-D Gaussian at a fixed sky position in each downloaded FITS image, measure the background RMS in a surrounding annulus, and save a results CSV. Rows where the fit fails are recorded with an empty flux column and can be plotted as upper limits.

The backend is chosen automatically at runtime (see [casda imfit — photometry backends](#casda-imfit--photometry-backends) below).

```bash
# single source — fit all FITS images under ./results/M31/images/
casda imfit -c "10.684 41.268" -n M31 -o ./results

# point at a custom image directory
casda imfit -c "10.684 41.268" -n M31 -o ./results --images-dir ./my_fits/

# batch from CSV
casda imfit --csv sources.csv -o ./results

# skip the automatic lightcurve plot
casda imfit -c "10.684 41.268" -n M31 -o ./results --no-plot

# tighter fit region, wider RMS annulus, μJy output
casda imfit -c "10.684 41.268" -n M31 -o ./results --radius 3 --rms-radius 45 --units uJy
```

Key options:

| Option | Default | Description |
|--------|---------|-------------|
| `--radius` | 5.0″ | Fit circle radius |
| `--rms-radius` | 30.0″ | Annulus outer radius for background RMS |
| `--scale` | 0.10 | Fractional flux-scale uncertainty |
| `--units` | mJy | Output flux units (`mJy` or `uJy`) |
| `--no-plot` | off | Skip the lightcurve PDF |

Results land in `OUTPUT/NAME/imfit/`.

#### casda imfit plot

Re-plot a lightcurve from an existing imfit results CSV without re-running the fits. Rows with a flux measurement are shown as scatter points (coloured by frequency when available); rows with no measurement are shown as downward-pointing triangles at 3 × RMS.

```bash
casda imfit plot --csv ./results/M31/imfit/M31_imfit_results.csv -n M31

# custom output directory and upper-limit threshold
casda imfit plot --csv ./results/M31/imfit/M31_imfit_results.csv -n M31 \
    -o ./figures --uplim-sigma 5
```

### casda cleanup

Remove large intermediate files (full-field FITS images and/or catalogue XMLs) while keeping all results (CSVs, cutouts, imfit outputs, plots).

```bash
# preview what would be removed
casda cleanup -n M31 -o ./results --dry-run

# remove images and catalogues for one source
casda cleanup -n M31 -o ./results

# keep images, only remove catalogues
casda cleanup -n M31 -o ./results --keep-images

# clean all sources in the output directory
casda cleanup -o ./results
```

## casda imfit — photometry backends

`casda imfit` fits a 2-D Gaussian to each downloaded FITS image at your source
position and measures the background RMS in a surrounding annulus.  It supports
two backends, selected automatically at runtime:

| Backend | When used | Platform |
|---------|-----------|----------|
| **CASA** (`casatasks`) | `casatasks` importable | Linux (PyPI wheels) |
| **astropy** (`astropy.modeling`) | fallback — always available | All (macOS, Linux, Windows) |

Both backends produce identical output columns and CSV schema.  The active
backend is printed at the start of every run so you always know which one is in
use.

**Installing the CASA backend (Linux only)**

```bash
pip install 'casda[casa]'
# or, in an editable install:
pip install -e '.[casa]'
```

CASA is not available as a pip package on macOS or Windows.  On those platforms
the astropy backend is used automatically — no extra steps needed.

**What the astropy backend does**

* Opens the FITS file with `astropy.io.fits`, squeezes degenerate Stokes/frequency axes.
* Converts your sky coordinates to pixel position via the image WCS.
* Fits a `Gaussian2D` model to a cutout around the source using
  `astropy.modeling.fitting.LevMarLSQFitter`.
* Measures background RMS as the standard deviation of pixels in an annulus
  (`2×radius` to `rms_radius` arcsec) around the source.
* Propagates fit, RMS, and flux-scale uncertainties in quadrature, identical to
  the CASA path.

## Updates:

#### October 14th 2026
The CLI (`cli-tools` branch) is now feature-complete and ready for wider testing.
- `casda lightcurves`, `casda images`, `casda imfit`, `casda cleanup`, and `casda credentials` are all fully documented above.
- `casda imfit` now works on macOS out of the box via an `astropy.modeling` fallback backend; CASA (`casatasks`) is an optional Linux-only extra (`pip install 'casda[casa]'`).
- `casda imfit plot` subcommand added: re-plot a lightcurve from an existing imfit results CSV without re-running fits. Non-detections (empty flux column) are shown as 3σ upper-limit arrows.
- `--no-plot` flag added to `casda imfit` (mirrors the existing flag on `casda lightcurves`).
- Parallel download threads increased to 8 by default across `casda lightcurves` and `casda images`.
- Credential injection fixed: OPAL credentials stored via `casda credentials save` are now correctly injected into astroquery, eliminating password prompts mid-run.

#### June 30th 2026
We have merged the working versions of the the `CASDA-Cutout` notebook and the `imfit_cli` script with the main branch.
These are still very much work in progress but we chose to make them accessible for beta testing.


#### June 26th 2026
Added a `imfit_cli.py` script to extract forced fits and limits:
- For use on image cutouts
- Currently available only on `dev` brahnch
- Will be implemented in the future to feed output into a combined lightcurve with the `Selavy` detections
  
#### March 20th 2026: 
We have further parallelised the downloads and retrieval and are in the process developing a CLI version and/or `uv` project to generate lightcurves. This will include integrating the image cutout tool for the purpose of forced "photometry" and limit measurements (currently in dev).

#### November 13th 2025: 
Added a `CASDA_Cutout` notebook to get image cutouts:
- Similarly optimised with batch staging and parallel downloads
- Currently available only on `dev` branch
- Single source processing has been tested
- Bulk downloads is untested **and not recommended** unless you want to fill up all of your storage

#### October 8th 2025: 
We overhauled `casda_download()` to optimise the process for speed.
- Chunked batch staging (10x faster than individual staging); `chunk_size` variable defined within function
- Parallel downloads for simultaneous file retrieval (3-5x faster); `max_workers` variable exposed to user
- Vectorized coordinate matching in `data_filter()` (2-5x faster)
- Dictionary lookups instead of DataFrame filtering (10-100x faster for lookups)
Overall: 20-100x speedup depending on your specific use case.
For benchmarking purposes you can set `optimized=False`.
