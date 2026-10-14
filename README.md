# CASDA
## a repository for CASDA-related tools

I've developed this notebook (based on code obtained from Dr. Laura Driessen, who got it from Dr. Minh Huynh) to download publicly available CASDA data around a given object's coordinates (or a .csv of objects/coordinates). [View Notebook](https://github.com/krose098/CASDA/blob/main/CASDA_Lightcurve_Public.ipynb)

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
