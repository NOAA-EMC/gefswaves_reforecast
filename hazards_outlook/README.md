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
cumulative distribution function. A predefined percentile, tuned to 87, is used to subsample the extreme tail of the
distribution, from which probabilities associated with given thresholds are computed. Following group meetings with
forecasters from OPC/NWS and considering the Beaufort scale and Saffir-Simpson Hurricane scale, the following levels
have been defined: wind speeds of 28, 34, 41, 48, 56, and 64 knots; significant wave heights of 4, 6, 9, and 14 meters;
and peak periods of 17, 20, and 22 seconds.
The probability maps have been generated on a daily basis since June, 2023, and the program has been run for
historical conditions related to significant case studies to validate the method across different meteorological events.
