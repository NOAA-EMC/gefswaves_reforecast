#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
extract_altimeter_fpos.py 
for fixed point positions

VERSION AND LAST UPDATE:
 v1.0  05/09/2022
 v1.1  04/10/2025

PURPOSE:
 Script to take altimeter tracks and collocate into fixed lat/lon positions.
 A total of 17 satellite missions are listed below. The period of each
  altimeter can be verified at:
  http://thredds.aodn.org.au/thredds/catalog/IMOS/SRS/Surface-Waves/Wave-Wind-Altimetry-DM00/catalog.html
  https://www.sciencedirect.com/science/article/pii/S0273117721000594
  https://ars.els-cdn.com/content/image/1-s2.0-S0273117721000594-gr1_lrg.jpg
 The AODN altimeter data can be downloaded with wfetchsatellite_AODN_Altimeter.sh from WW3-tools:
  https://github.com/NOAA-EMC/WW3-tools

USAGE:
 This program processes Altimeter data from AODN and collocates it into a specified
  ID/Lat/Lon positions, which is provided as an input txt file.
 Altimeters must have been downloaded beforehand (see wfetchsatellite_AODN_Altimeter.sh)
 Check the pre-selected parameters below for the altimeter collocation
 Input arguments are:
   1) path/name.txt of point positions (ID Lat Lon)
   2) Initial date (YYYYMMDD)
   3) Final date (YYYYMMDD)
   4) Delta time (temporal resolution), in hours
 Example (from linux terminal command line):
   nohup python3 extract_altimeter_fpos.py points_sel.txt 20201001 20250101 6 >> nohup_extract_Altimeter.out 2>&1 &

OUTPUT:
 Netcdf file containing the collocated altimeter data into fixed lat/lons.
 Hs: significant wave height
 U10: 10-meter wind speed

DEPENDENCIES:
 See setup.py and the imports below.
 AODN altimeter data previously downloaded (see wfetchsatellite_AODN_Altimeter.sh)

AUTHOR and DATE:
 05/09/2022: Ricardo M. Campos, first version.
 04/10/2025: Ricardo M. Campos, Using fixed position given by a text file. Save output in netcdf format.

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib.mlab import *
# from pylab import *
import statistics
import netCDF4 as nc
import pandas as pd
import time
import timeit
from time import strptime
from calendar import timegm
import datetime as dt
import sys
import warnings; warnings.filterwarnings("ignore")
# netcdf format
fnetcdf="NETCDF4"


# Distance (km) between two points. A function that will be used during collocation loop.
def distance(lat1, lon1, lat2, lon2):
    # convert decimal degrees to radians
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    # haversine formula
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = math.sin(dlat/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2)**2
    c = 2 * math.asin(math.sqrt(a))
    r = 6371  # radius of Earth in kilometers
    return c * r

if __name__ == "__main__":

    # Input text file name (with path if not in the same local dir)
    fname=str(sys.argv[1]) # point_name(fixed) lat lon to be extracted.
    # fname='/work/noaa/marine/ricardo.campos/work/analysis/3assessments/fuzzy_verification/points_sel.txt'
    idate=str(sys.argv[2]); fdate=str(sys.argv[3]) # YYYYMMDD
    dt=int(sys.argv[4])

    start = timeit.default_timer()

    # Time array
    ftime = np.arange(timegm( strptime(idate+"0000",'%Y%m%d%H%M')),timegm( strptime(fdate+"0000",'%Y%m%d%H%M'))+1.,dt*3600.).astype('double')

    # Directory where AODN altimeter data is saved, downloaded using wfetchsatellite_AODN_Altimeter.sh
    dirs='/work/noaa/marine/ricardo.campos/data/AODN/altimeter'
    # Minimum distance (km) from the coast
    mindfc=30. # in Km
    # Maximum distance (m) for the collocation average
    dlim=37000. # in m
    # Maximum temporal distance (s) for the collocation average
    maxti=3600.
    # power of initial array 10**pia (size) that will be used to allocate satellite data (faster than append)
    pia=10
    # Satellite missions available at AODN dataset.
    sdname=np.array(['JASON3','JASON2','CRYOSAT2','JASON1','HY2','HY2B','SARAL','SENTINEL3A','ENVISAT','ERS1','ERS2','GEOSAT','GFO','TOPEX','SENTINEL3B','CFOSAT','SENTINEL6A'])
    sname=np.array(['JASON-3','JASON-2','CRYOSAT-2','JASON-1','HY-2','HY-2B','SARAL','SENTINEL-3A','ENVISAT','ERS-1','ERS-2','GEOSAT','GFO','TOPEX','SENTINEL-3B','CFOSAT','SENTINEL-6A'])

    # Altimeter Quality Control parameters
    max_swh_rms = 1.5  # Max RMS of the band significant wave height
    max_sig0_rms = 0.8 # Max RMS of the backscatter coefficient
    max_swh_qc = 2.0 # Max SWH Ku band quality control
    hsmax=20.; wspmax=80.
    min_swh_numval = np.array([17,17,17,17,17,17,17,17,17,17,17,-np.inf,3,7,17,-np.inf,17])

    # --- INPUT ARRAY Information ---

    # read points of interest
    # ds_points = pd.read_csv(fname,comment='#',delimiter=r"\s+")
    ds_points = pd.read_csv(fname,comment='#',delimiter=" ")
    bid=np.array(ds_points['ID'].values).astype('str'); blat=np.array(ds_points['Lat'][:]); blon=np.array(ds_points['Lon'][:]); glon=np.copy(blon); glon[glon<0]=glon[glon<0]+360.

    # initialize arrays
    fhskcal=np.zeros((len(bid),ftime.shape[0]),'f')*np.nan; fwndcal=np.zeros((len(bid),ftime.shape[0]),'f')*np.nan

    for i in range(0,len(bid)):

        print(" ---- Point "+bid[i]+" ---- ")
        # Sat files (squares) considering the lat lon of interest, for the AODN file names
        auxlat=np.array(np.arange(np.floor(np.nanmin(blat[i]))-1.,np.ceil(np.nanmax(blat[i]))+1.,1)).astype('int')
        auxlon=np.array(np.arange(np.floor(np.nanmin(glon[i]))-1.,np.ceil(np.nanmax(glon[i]))+1.,1)).astype('int')

        # Read and allocate satellite data into arrays
        ast=np.double(np.zeros((10**pia),'d')); adfc=np.zeros((10**pia),'f')
        aslat=np.zeros((10**pia),'f'); aslon=np.zeros((10**pia),'f')
        ahskcal=np.zeros((10**pia),'f'); awndcal=np.zeros((10**pia),'f')
        asig0knstd=np.zeros((10**pia),'f'); aswhknobs=np.zeros((10**pia),'f')
        aswhknstd=np.zeros((10**pia),'f'); aswhkqc=np.zeros((10**pia),'f')
        ii=0
        for s in range(0,np.size(sdname)):
            for j in auxlat:
                for k in auxlon:

                    if j>=0:
                        hem='N'
                    else:
                        hem='S'

                    try: 
                        fu=nc.Dataset(dirs+'/'+sdname[s]+'/IMOS_SRS-Surface-Waves_MW_'+sname[s]+'_FV02_'+str(np.abs(j)).zfill(3)+hem+'-'+str(k).zfill(3)+'E-DM00.nc')
                    except:
                        print(dirs+'/'+sdname[s]+'/IMOS_SRS-Surface-Waves_MW_'+sname[s]+'_FV02_'+str(np.abs(j)).zfill(3)+hem+'-'+str(k).zfill(3)+'E-DM00.nc does not exist'); vai=0
                    else:
                        st=np.double(fu.variables['TIME'][:]*24.*3600.+float(timegm( time.strptime('1985010100', '%Y%m%d%H') )))
                        if (np.size(st)>10) and (np.nanmax(st)>np.nanmin(ftime)) and (np.nanmin(st)<np.nanmax(ftime)):
                            slat=fu.variables['LATITUDE'][:]
                            slon=fu.variables['LONGITUDE'][:]
                            wndcal=fu.variables['WSPD_CAL'][:]
                            try: 
                                hskcal=fu.variables['SWH_KU_CAL'][:]
                                sig0knstd=fu.variables['SIG0_KU_std_dev'][:]
                                swhknobs=fu.variables['SWH_KU_num_obs'][:]
                                swhknstd=fu.variables['SWH_KU_std_dev'][:]
                                swhkqc=fu.variables['SWH_KU_quality_control'][:]
                                dfc=fu.variables['DIST2COAST'][:]
                            except:
                                print(' error reading KU, pick KA')
                                hskcal=fu.variables['SWH_KA_CAL'][:]
                                sig0knstd=fu.variables['SIG0_KA_std_dev'][:]
                                swhknobs=fu.variables['SWH_KA_num_obs'][:]
                                swhknstd=fu.variables['SWH_KA_std_dev'][:]
                                swhkqc=fu.variables['SWH_KA_quality_control'][:]
                                dfc=swhkqc*0.+999.

                            if ii+np.size(st) <= ast.shape[0] :
                                if (st.shape[0]==wndcal.shape[0]) and (slat.shape[0]==slon.shape[0]) and (hskcal.shape[0]==wndcal.shape[0]) :    
                                    ast[ii:ii+st.shape[0]]=np.array(st).astype('double')
                                    aslat[ii:ii+st.shape[0]]=np.array(slat).astype('float')
                                    aslon[ii:ii+st.shape[0]]=np.array(slon).astype('float')
                                    ahskcal[ii:ii+st.shape[0]]=np.array(hskcal).astype('float')
                                    awndcal[ii:ii+st.shape[0]]=np.array(wndcal).astype('float')
                                    asig0knstd[ii:ii+st.shape[0]]=np.array(sig0knstd).astype('float')
                                    aswhknobs[ii:ii+st.shape[0]]=np.array(swhknobs).astype('float')
                                    aswhknstd[ii:ii+st.shape[0]]=np.array(swhknstd).astype('float')
                                    aswhkqc[ii:ii+st.shape[0]]=np.array(swhkqc).astype('float')
                                    adfc[ii:ii+st.shape[0]]=np.array(dfc).astype('float')
                                    ii=ii+st.shape[0]

                            else:
                                sys.exit('Small array to allocate the satellite data! Increase the power of initial array (pia)')

                            del st,slat,slon,hskcal,wndcal,sig0knstd,swhknobs,swhknstd,swhkqc,dfc
                            fu.close(); del fu

            print('    done reading and allocating satellite data '+sdname[s])

        del ii
        print('--- Done reading and allocating satellite data ---')
        print('  Running collocation ... ')
        # Simplified Quality Control Check ----
        adatemin=np.nanmin(ftime)-3600.; adatemax=np.nanmax(ftime)+3600.
        indq = np.where( (adfc>=mindfc) & (aswhknstd<=max_swh_rms) & (asig0knstd<=max_sig0_rms) & (aswhknobs>=min_swh_numval[s]) & (aswhkqc<=max_swh_qc) & (ahskcal>0.1) & (ahskcal<hsmax) & (awndcal>0.2) & (awndcal<wspmax) & (ast>=adatemin) & (ast<=adatemax) )     
        del asig0knstd,aswhknobs,aswhknstd,aswhkqc,adatemin,adatemax,adfc

        if np.size(indq)>10:
            ii=0 # apply quality control
            ast=np.double(np.copy(ast[indq[0]]))
            aslat=np.copy(aslat[indq[0]]); aslon=np.copy(aslon[indq[0]])
            aslon[aslon>180.]=aslon[aslon>180.]-360.
            ahskcal=np.copy(ahskcal[indq[0]]); awndcal=np.copy(awndcal[indq[0]])

            for t in range(0,ftime.shape[0]):
                indt = np.where( abs(ast[:]-ftime[t]) < maxti )
                if np.size(indt)>1:
                    # print(repr(t))

                    # distance
                    sdist=np.zeros((np.size(indt)),'f')*np.nan
                    for j in range(0,np.size(indt)):
                        sdist[j]=distance(aslat[indt][j],aslon[indt][j],blat[i],blon[i])

                    # within search radius
                    inds=np.where((sdist<=dlim/1000.)&(ahskcal[indt]<hsmax))
                    if np.size(inds)>1:
                        fhskcal[i,t]=np.float(np.nanmean(ahskcal[indt][inds]))
                        fwndcal[i,t]=np.float(np.nanmean(awndcal[indt][inds]))

                    del inds,sdist

                del indt
                # print('  -- collocation '+repr(t)+' . '+repr(ftime.shape[0]))

            del ast, aslat, aslon, ahskcal, awndcal, indq

            # Final quality control (double check)
            ind=np.where( (fwndcal[i,:]<0.01) | (fwndcal[i,:]>wspmax) )
            if np.any(ind):
                fwndcal[i,:][ind[0]]=np.nan; del ind

            ind=np.where( (fhskcal[i,:]<0.1) | (fhskcal[i,:]>hsmax) )
            if np.any(ind):
                fhskcal[i,:][ind[0]]=np.nan; del ind

            print(" Total altimeter collocated records (point "+bid[i]+") : "+str(np.size(ind))); print(' ')

        else:
            print(" No satellite records within the given time/date range and quality control parameters selected. Point "+bid[i])

        print(" ---- ")

    print(' ')
    # Save netcdf output file 
    ncfile = nc.Dataset("Altimeter.PointExtract_"+idate+"to"+fdate+".nc", "w", format=fnetcdf)
    ncfile.history="AODN Altimeter data extracted for fixed point outputs."
    # create  dimensions
    ncfile.createDimension('points', bid.shape[0] )
    ncfile.createDimension('time', len(ftime) )
    # create variables.
    vt = ncfile.createVariable('time',np.dtype('float64').char,('time'))
    vbid = ncfile.createVariable('ID',np.dtype('a25'),('points'))
    vblat = ncfile.createVariable('lat',np.dtype('float32').char,('points'))
    vblon = ncfile.createVariable('lon',np.dtype('float32').char,('points'))
    # results
    vu10 = ncfile.createVariable('u10',np.dtype('float32').char,('points','time'))
    vhs = ncfile.createVariable('hs',np.dtype('float32').char,('points','time'))
    # Assign units
    vt.units = 'seconds since 1970-01-01T00:00:00+00:00'
    vblat.units = 'degrees_north' ; vblon.units = 'degrees_east'
    vu10.units = 'm/s'; vhs.units = 'm'
    # Allocate Data
    vt[:]=ftime[:]; vbid[:]=bid[:]
    vblat[:]=np.array(blat[:]).astype('float'); vblon[:]=np.array(blon[:]).astype('float')
    vu10[:,:]=fwndcal[:,:]; vhs[:,:]=fhskcal[:,:]
    ncfile.close()
    print("Done. Netcdf ok. New file saved: Altimete.PointExtract_"+idate+"to"+fdate+".nc")

    stop = timeit.default_timer()
    print('Concluded in '+repr(int(round(stop - start,0)))+' seconds')
    print(" ------------------ "); print(" ")


