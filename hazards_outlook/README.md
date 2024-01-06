About the ensemble wave forecast GEFSv12:
The wave component from GEFS is generated using the numerical wave model WAVEWATCH III version 7.0.
The wave forecasts are on a grid of spatial resolution of 0.25° x 0.25° and products are available on a 
temporal resolution of 3 hours for the first 10 days and 6 hours up to 16 days. 
The ensemble consists of 30 perturbed members plus the control member. See:
https://noaa-nws-gefswaves-reforecast-pds.s3.amazonaws.com/index.html
https://github.com/NOAA-EMC/WW3-tools/tree/develop/opforecast/get_gefsWaves.sh
https://github.com/NOAA-EMC/gefswaves_reforecast/tree/main/download_and_plot

The algorithm for probability maps begins by pooling forecast data into 2° x 2° bins, centered at each grid point,
over a 7-day time window associated with the second week of the forecast. This process is repeated for the entire
ensemble member set, which is then combined and reshaped into a single large array, forming an empirical
cumulative distribution function. A predefined percentile is used to subsample the extreme tail of the
distribution, from which probabilities associated with given thresholds are computed. Considering the Beaufort scale 
and Saffir-Simpson Hurricane scale, the following levels have been defined: 
wind speeds of 41, 48, 56, and 64 knots; significant wave heights of 4, 6, 9, and 14 meters
The probability maps have been generated on a daily basis since June, 2023, and the program has been run for
historical conditions related to significant case studies to validate the method across different meteorological events.

probmaps_gefs.yaml:
It is a configuration file containing the paths, to be edited, and technical information for
the probability maps.

download_GEFSwaves.sh:
Script to download the NOAA Global Ensemble Forecast System (GEFS), Wave 
Forecast. It downloads the field outputs in grib2 format.

probmaps_gefs.sh:
Script to check if GEFS files are available, and then runs the python code
to generate the probabilies. This is the script that should be in the crontab,
running daily.

probmaps_gefs.py:
Script with the statistical analysis and plots, generating the probability maps.

