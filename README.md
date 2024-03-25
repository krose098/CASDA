# CASDA
repository for CASDA-related tools

I've developed this notebook (based on code obtained from Dr. Laura Driessen, who got it from Dr. Minh Huynh) to download publicly available CASDA data around a given object's coordinates (or a .csv of objects/coordinates).

All you need to do is have an OPAL login, add your source details, (pray that CASDA is working), and it will:
1) Download the component catalogues into a newly created, named folder (you can play around with the code and download other data products too)
2) Match components to within a given sky separation (default 5")
3) Produce a lightcurve from all public ASKAP data
4) Save the lightcurve data points to a .csv.

A couple of small things to keep in mind with this data:
* The matches are to Selavy component source catalogues (within 5") of the object coordinates and Selavy has a 5sigma detection threshold.
* There was a recent bug fix due to a change in the CASDA file naming convention. It is possible that some files are missed by my search if they are unconventionally named.
* In some cases CASDA may store multiple versions of a single observation -- especially from pilot programs -- if the data was reprocessed. This may produce, for example, two detections on the same day with slightly different fluxes.
* The columns in each .csv should be self-explanatory (feel free to reach out if not) with the exception of flux_err_quad which is a more conservation flux uncertainty I defined for the lightcurves as: ((flux_err^2) + (rms^2) + (0.06*flux)^2)^0.5
