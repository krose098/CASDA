# CASDA
## a repository for CASDA-related tools

I originally developed this notebook (based on code obtained from Dr. Laura Driessen, who got it from Dr. Minh Huynh) to download publicly available CASDA data around a given object's coordinates (or a .csv of objects/coordinates). [View Notebook](https://github.com/krose098/CASDA/blob/main/CASDA_Lightcurve_Public.ipynb). See the `cli-tools` for installation and usage guide for the command line interface (CLI) version.


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
- Parallel downloads for simultaneous file retrieval (3-5x faster); `max_workers` variable exposed to user; this should be tested with low numbers (2-4) prior to accelerating with (8-12) workers, depending on your system.
- Vectorized coordinate matching in `data_filter()` (2-5x faster)
- Dictionary lookups instead of DataFrame filtering (10-100x faster for lookups)
Overall: 20-100x speedup depending on your specific use case.
For benchmarking purposes you can set `optimized=False`.
