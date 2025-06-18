#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
extractPoint_makeHindcast_GEFSreforecast.py

VERSION AND LAST UPDATE:
 v1.0  08/26/2023

PURPOSE:
 Script to extract point outputs from ww3 field outputs (wave ensemble, GEFS)
  and build a continuous hindcast data (append consecutive 24h slices).
 Time-series of several ww3 variables (lvars) are produced for the given points.
 Each run extracts one month only, from GEFS reforecast data:
  download_AWSreforecast_Field.sh
  https://noaa-nws-gefswaves-reforecast-pds.s3.amazonaws.com/index.html
 A list of lon lat buoyID must be given, in text format.

USAGE:
 This program reads grib2 field outputs from GEFSv12 reforecast (all ensemble members)
  and extract information for a given lon/lat position.
 The reforecast must have been previously downloaded. See download_AWSreforecast_Field.sh
 The variables outpath, gefs_path, buoypfile below need to be edited.
 Two input arguments: Year Month
 Example (from linux terminal command line):
  nohup python3 extractPoint_makeHindcast_GEFSreforecast.py 2010 1 >> nohup_extractPoint_makeHindcast.out 2>&1 &
 Or you can submit a job to run it (Orion supercomputer), using jextractPoint_makeHindcast.sh
 Multiple jobs can be submitted using jobmanager.sh

OUTPUT:
 Netcdf file GEFS.PointExtract.Day1.*.nc containing the time-series of 
  several variables (lvars) for the given point positions, covering a 
  period of one month.

DEPENDENCIES:
 See setup.py and the imports below.
 GEFSv12 reforecast data (see download_AWSreforecast_Field.sh)

AUTHOR and DATE:
 08/26/2023: Ricardo M. Campos, first version.

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

import matplotlib
matplotlib.use('Agg')
from matplotlib.mlab import *
import numpy as np
import sys
import pandas as pd
import xarray as xr
import netCDF4 as nc
import time
from time import strptime
from calendar import timegm
from datetime import datetime, timedelta
import warnings; warnings.filterwarnings("ignore")
# netcdf format
fnetcdf="NETCDF4"

# fixed variables
# output path where the final netcdf file will be saved.
outpath="/work/noaa/marine/ricardo.campos/work/analysis/1preproc/extrai_point_from_grib/"
# model path where the ww3 field output is saved.
gefs_path="/work/noaa/marine/ricardo.campos/work/analysis/GEFSv12/"
# ensemble members
ensm=['c00','p01','p02','p03','p04']
# variable names (exactly as written in the ww3 field output)
lvars=['swh','u','v','mwp','perpw','mpww','mwd','dirpw','shww','wvdir']
# position (lon lat) and buoy ID, saved in a text file. Longitude from -180 to 180.
buoypfile="/work/noaa/marine/ricardo.campos/work/analysis/1preproc/extrai_point_from_grib/bstations_GEFSv12WW3grid.txt"
# Example of one line
# -160.794  53.969  46075

# time resolution in hours of the output time array
tres=3.

# Input arguments
year=int(sys.argv[1])
month=int(sys.argv[2])

# Calculate the next month
next_month_date = datetime(int(year),int(month), 1) + timedelta(days=31)
# time resolution
indl=int(24/tres)
# Time intervall and array of reference
firstdate=repr(year)+str(month).zfill(2)+'01'; lastdate=repr(next_month_date.year)+str(next_month_date.month).zfill(2)+'01'
# time array
ftime=np.array(np.arange(timegm( strptime(firstdate, '%Y%m%d') ),timegm( strptime(lastdate, '%Y%m%d') )+1,tres*3600.)).astype('double')
ftime=ftime[1::]
# cycle array. One cycle per day. Edit here if you have more cycles per day.
ctime=np.array(np.arange(timegm( strptime(firstdate, '%Y%m%d') ),timegm( strptime(lastdate, '%Y%m%d') ),24*3600.)).astype('double')

# READ BUOY POSITIONS
bpos = pd.read_csv(buoypfile, delim_whitespace=True)
latb=np.array(bpos['lat']); lonb=np.array(bpos['lon']); bid=np.array(bpos['buoyID'])
lonb[lonb>180]=lonb[lonb>180]-360.
# -----------

winterp=np.zeros((len(lvars),len(bid),len(ensm),len(ftime)),'f')*np.nan
c=0
for j in range(0,len(ensm)):
    t=0
    for i in range(0,len(ctime)):
        fname=gefs_path+"gefs.wave."+repr(time.gmtime(ctime[i])[0])+str(time.gmtime(ctime[i])[1]).zfill(2)+str(time.gmtime(ctime[i])[2]).zfill(2)+"."+ensm[j]+".global.0p25.grib2"
        try:
            f = xr.open_dataset(fname, engine='cfgrib')

            if c==0:
                # Origin
                latm=np.array(f['latitude'].values); lonm=np.array(f['longitude'].values)
                lonm[lonm>180]=lonm[lonm>180]-360.
                indlat=[]; indlon=[]
                for k in range(0,len(bid)):
                    indlat=np.append(indlat,np.min(np.where(abs(latm-latb[k])==np.min(abs(latm-latb[k])))[0]))
                    indlon=np.append(indlon,np.min(np.where(abs(lonm-lonb[k])==np.min(abs(lonm-lonb[k])))[0]))

                indlat=np.array(indlat).astype('int'); indlon=np.array(indlon).astype('int')

            # extract points
            for v in range(0,len(lvars)):
                winterp[v,:,j,t:t+indl] = np.array(f[lvars[v]].isel(step=np.arange(0,indl)).values[:,indlat,indlon]).T

        except:
            print(" Error: Cannot open and read "+fname)

        else:
            t=t+indl; c=c+1
            print(fname)
            f.close(); del f, fname


# Save netcdf output file
ncfile = nc.Dataset(outpath+"GEFS.PointExtract.Day1."+firstdate+"to"+lastdate+".nc", "w", format=fnetcdf)
ncfile.history="Wave parameters extracted for buoy points, using bstations_GEFSv12WW3grid.txt. Forecast lead time 0-24h."
# create  dimensions. 2 Dimensions
ncfile.createDimension('variabels',len(lvars))
ncfile.createDimension('buoy_points', len(bid))
ncfile.createDimension('ensemble_member', len(ensm))
ncfile.createDimension('time', len(ftime) )
# create variables.
vt = ncfile.createVariable('time',np.dtype('float64').char,('time'))
vbid = ncfile.createVariable('buoyID',np.dtype('a25'),('buoy_points'))	
vlat = ncfile.createVariable('latitude',np.dtype('float32').char,('buoy_points'))
vlon = ncfile.createVariable('longitude',np.dtype('float32').char,('buoy_points'))
vensm = ncfile.createVariable('ensemble_member',np.dtype('a25'),('ensemble_member'))
vlvars = ncfile.createVariable('variables_names',np.dtype('a25'),('variabels'))
# results
vwinterp = ncfile.createVariable('gefs_ww3',np.dtype('float32').char,('variabels','buoy_points','ensemble_member','time'))
# Assign units
vt.units = 'seconds since 1970-01-01T00:00:00+00:00'
# Allocate Data
vt[:]=ftime[:]; vbid[:]=bid[:]; vlat[:]=latb[:]; vlon[:]=lonb[:]
vensm[:]=np.array(ensm[:]); vlvars[:]=np.array(lvars[:])
vwinterp[:,:,:,:]=winterp[:,:,:,:]
ncfile.close()
print(' ')
print('Done. Netcdf ok. New file saved: GEFS.PointExtract.'+firstdate+'to'+lastdate+'.nc')

