#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
buildfuzzydataset.py

VERSION AND LAST UPDATE:
 v1.0  09/10/2023
 v1.1  05/09/2024

PURPOSE:
 Python script to build a dataset with GEFS forecast, GEFS forecast (hindcasted 24-h slices),
  and NDBC buoy data. It uses the NDBC buoy point positions to select the nearest model point and
  neighbouring points using a window (wl size) centered at the NDBC buoy point.
  This is intended for fuzzy verification.

USAGE:
 This code is run for one forecast cycle only. A shell script can execute this script multiple times.
 python3 buildfuzzydataset.py 20220701

OUTPUT:
 Netcdf file GEFS.BUOY.PointExtract.YYYYMMDD.nc saved in the given outpath.

DEPENDENCIES:
 See setup.py and the imports below.

AUTHOR and DATE:
 09/10/2023: Ricardo M. Campos, first version.
 05/09/2024: Ricardo M. Campos. GDAS removed. Now reading post-processed netcdf file instead of .grib2 (much faster).

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

import numpy as np
import sys
import pandas as pd
import xarray as xr
import netCDF4 as nc
import time
import timeit
from time import strptime
from calendar import timegm
import wread
import warnings; warnings.filterwarnings("ignore")
# netcdf format
fnetcdf="NETCDF4"


def read_gefs(date,mvar,path_gefs):

    # READ WW3 Ensemble Forecast files. Appending forecast days (ensemble members).
    c=0
    for enm in range(0,nenm):

        fname=path_gefs+"gefs.wave."+date+"."+str(enm).zfill(2)+".global.0p25.nc"
        if c==0:
            ds = xr.open_dataset(fname)
            wdate = pd.to_datetime(np.atleast_1d(np.array(ds.time.values)))
            wtime = np.atleast_1d(np.array(ds.time.values))

            lat = ds['latitude'].values[:]; lon = ds['longitude'].values[:]

            fmod = np.zeros((len(mvar),wtime.shape[0],nenm,lat.shape[0],lon.shape[0]),'f')*np.nan
            ds.close(); del ds

        f=nc.Dataset(fname)
        for i in range(0,len(mvar)):
            fmod[i,:,enm,:,:] = np.array(f.variables[mvar[i]][:,:,:])

        f.close(); del f
        print("    read_gefs enm "+repr(enm)+" ok")
        c=c+1

    wtime = np.array( (wtime.astype('int64') // 10**9) ).astype('double')

    return np.array(fmod), np.array(lat), np.array(lon), wtime, wdate


def build_gefs_hindcast(date,mvar,path_gefs):

    # cycle time
    ctime=np.array(np.arange(0,384+1,24)*3600.+timegm( strptime(date,'%Y%m%d') )).astype('double')
    auxltime=np.array(np.arange(0,24,6)).astype('int'); auxltimef=np.array([0]).astype('int')

    c=0; tt=0
    for ct in range(0,len(ctime)):

        if ct==int(len(ctime)-1):
            auxlt=auxltimef
        else:
            auxlt=auxltime

        for enm in range(0,nenm):
            cdate=str(time.gmtime(ctime[ct])[0])+str(time.gmtime(ctime[ct])[1]).zfill(2)+str(time.gmtime(ctime[ct])[2]).zfill(2)
            fname=path_gefs+"gefs.wave."+cdate+"."+str(int(enm)).zfill(2)+".global.0p25.nc"
            # print(" build_gefs_hindcast "+fname)

            if c==0:
                ds = xr.open_dataset(fname)
                wdate = pd.to_datetime(np.atleast_1d(np.array(ds.time.values)))
                wtime = np.atleast_1d(np.array(ds.time.values))
                lat = ds['latitude'].values[:]; lon = ds['longitude'].values[:]
                fmod=np.zeros((len(mvar),len(wtime),nenm,len(lat),len(lon)),'f')*np.nan
                ds.close(); del ds

            f=nc.Dataset(fname)
            for i in range(0,len(mvar)):
                fmod[i,tt:tt+len(auxlt),enm,:,:] = np.array(f.variables[mvar[i]][0:len(auxlt),:,:])

            f.close(); del f
            c=c+1

        print("    build_gefs_hindcast cycle "+repr(ct)+" ok")
        tt=tt+len(auxlt)

    wtime = np.array( (wtime.astype('int64') // 10**9) ).astype('double')

    return fmod, np.array(lat), np.array(lon), wtime, wdate


def read_buoy(date,bid,anh,blat,blon,path_buoy):

    auxltime=np.array(np.arange(0,384+1,6)).astype('int')
    atime=np.array(auxltime*3600.+timegm( strptime(date,'%Y%m%d') )).astype('double')
    # year
    datestrg1=str(time.gmtime(atime[0])[0]); datestrg2=str(time.gmtime(atime[-1])[0])
    pskip = 0
    if datestrg1 == datestrg2 :

        try: 
            bdata=wread.tseriesnc_ndbc(path_buoy+str(bid)+"h"+datestrg1+".nc",anh=anh)
            bt=np.array(bdata['time']).astype('double')
            fbdata=np.array([bdata['wind_spd']])
            fbdata=np.append(fbdata,[bdata['hs']],axis=0)
            fbdata=np.append(fbdata,[bdata['tp']],axis=0)
            ndbclat=float(bdata['latitude']); ndbclon=float(bdata['longitude'])
            del bdata
        except:
            print("   No file "+path_buoy+str(bid)+"h"+datestrg1+".nc")
            fdata = np.zeros((3,len(atime)),'f')*np.nan
            ndbclat=blat; ndbclon=blon
            pskip = 1
        else:
            print("    read_buoy "+str(bid)+" "+str(date)+" ok")

    else:

        try: 
            bdata1=wread.tseriesnc_ndbc(path_buoy+str(bid)+"h"+datestrg1+".nc",anh=anh)
            bdata2=wread.tseriesnc_ndbc(path_buoy+str(bid)+"h"+datestrg2+".nc",anh=anh)
            bt=np.array(np.append(bdata1['time'],bdata2['time'])).astype('double')
            fbdata=np.array([np.append(bdata1['wind_spd'],bdata2['wind_spd'])])
            fbdata=np.append(fbdata,[np.append(bdata1['hs'],bdata2['hs'])],axis=0)
            fbdata=np.append(fbdata,[np.append(bdata1['tp'],bdata2['tp'])],axis=0)
            ndbclat=float(bdata1['latitude']); ndbclon=float(bdata1['longitude'])
            del bdata1,bdata2
        except:
            print("   No file "+path_buoy+str(bid)+"h"+datestrg1+".nc and/or "+path_buoy+str(bid)+"h"+datestrg2+".nc")
            fdata = np.zeros((3,len(atime)),'f')*np.nan
            ndbclat=blat; ndbclon=blon
            pskip = 1
        else:
            print("    read_buoy "+str(bid)+" "+str(date)+" ok")

    if pskip == 0:
        fdata=np.zeros((3,len(atime)),'f')*np.nan
        for i in range(0,len(atime)):
            indt=np.where(np.abs(atime[i]-bt)<3600)
            if np.size(indt)>0:
                indt=indt[0]
                for j in range(0,fbdata.shape[0]):
                    fdata[j,i]=np.nanmean(fbdata[j,indt])

                del indt

    return fdata, ndbclat, ndbclon, atime
    del pskip


if __name__ == "__main__":

    # Input argument (cycle time only):
    date=str(sys.argv[1])
    print(" "); print("    ------- "+date+" ------- ")

    # start time
    start = timeit.default_timer()

    # max spatial window radius(diameter/2)
    wl=3 # (°)
    # number of ensemble members
    nenm=int(30+1)
    # Data paths
    path_gefs="/work/noaa/marine/ricardo.campos/data/archiveOPruns/GEFSv12Waves_AWS/netcdf/"
    path_buoy="/work/noaa/marine/ricardo.campos/data/buoys/NDBC/ncformat/wparam/"
    outpath="/work/noaa/marine/ricardo.campos/work/analysis/3assessments/fuzzy_verification/data/"
    # Read Buoy points
    ds_buoy = pd.read_csv('/work/noaa/marine/ricardo.campos/work/analysis/3assessments/fuzzy_verification/buoys_sel.txt',comment='#',delimiter=r"\s+")

    # variable names
    mvar=np.array(["WIND_surface","HTSGW_surface","PERPW_surface"])

    # GEFS forecast
    gefs_forecast, gflat, gflon, gftime, gfdate = read_gefs(date,mvar,path_gefs)
    # GEFS sliced, build hindcast
    gefs_hindcast, ghlat, ghlon, ghtime, ghdate = build_gefs_hindcast(date,mvar,path_gefs)

    for i in range(0,len(ds_buoy)):
        bid=ds_buoy['buoyID'][i]; anh=ds_buoy['AnH'][i]
        blat=float(ds_buoy['Lat'][i]); blon=float(ds_buoy['Lon'][i])

        # buoy lon standard
        if blon<0:
            blon=blon+360

        # Buoy data
        buoy_data, bdlat, bdlon, bdtime = read_buoy(date,bid,anh,blat,blon,path_buoy)

        dlat=int(np.round(wl/np.mean(np.diff(gflat))))
        indlat=np.min(np.where(np.abs(gflat-blat) == np.min(np.abs(gflat-blat))))
        indlat=np.array(np.arange(int(indlat-dlat),int(indlat+dlat+1))).astype('int')

        dlon=int(np.round(wl/np.mean(np.diff(gflon))))
        indlon=np.min(np.where(np.abs(gflon-blon) == np.min(np.abs(gflon-blon))))
        indlon=np.array(np.arange(int(indlon-dlon),int(indlon+dlon+1))).astype('int')

        if (np.array_equal(gflat,ghlat)==True) and (np.array_equal(gflon,ghlon)==True) and (np.array_equal(gftime,ghtime)==True):
            if i==0:

                fgefs_forecast = np.array([gefs_forecast[:,:,:,indlat,:][:,:,:,:,indlon]])
                fgefs_hindcast = np.array([gefs_hindcast[:,:,:,indlat,:][:,:,:,:,indlon]])
                fbuoy_data = np.array([buoy_data])
                fgftime = np.array(gftime)
                fgflat = np.array([gflat[indlat]])
                fgflon = np.array([gflon[indlon]])

            else:

                fgefs_forecast = np.append(fgefs_forecast,np.array([gefs_forecast[:,:,:,indlat,:][:,:,:,:,indlon]]),axis=0)
                fgefs_hindcast = np.append(fgefs_hindcast,np.array([gefs_hindcast[:,:,:,indlat,:][:,:,:,:,indlon]]),axis=0)
                fbuoy_data = np.append(fbuoy_data,np.array([buoy_data]),axis=0)
                fgflat = np.append(fgflat,np.array([gflat[indlat]]),axis=0)
                fgflon = np.append(fgflon,np.array([gflon[indlon]]),axis=0)

        del buoy_data,dlat,dlon,indlat,indlon
        print(" "); print(" --- Point "+repr(bid)+" Ok."); print(" ")

    # keep model longitude with the same standard as the NDBC lon (-180 to 180)
    fbmlon=np.array(ds_buoy['Lon']).astype('float'); fbmlon[fbmlon>180]=fbmlon[fbmlon>180]-360.

    # Save netcdf output file
    ncfile = nc.Dataset(outpath+"GEFS.BUOY.PointExtract."+date+".nc", "w", format=fnetcdf)
    ncfile.history="Data extracted from GEFS, for NDBC buoy positions and neighbouring points."
    # create  dimensions.
    ncfile.createDimension('buoy_points', len(ds_buoy))
    ncfile.createDimension('ensemble_member', nenm)
    ncfile.createDimension('time', len(fgftime) )
    ncfile.createDimension('lat', fgflat.shape[1] )
    ncfile.createDimension('lon', fgflon.shape[1] )
    # create variables.
    vt = ncfile.createVariable('time',np.dtype('float64').char,('time'))
    vbid = ncfile.createVariable('buoyID',np.dtype('a25'),('buoy_points'))
    vblat = ncfile.createVariable('buoy_lat',np.dtype('float32').char,('buoy_points'))
    vblon = ncfile.createVariable('buoy_lon',np.dtype('float32').char,('buoy_points'))
    vensm = ncfile.createVariable('ensemble_member',np.dtype('int').char,('ensemble_member'))
    vlat = ncfile.createVariable('lat',np.dtype('float32').char,('buoy_points','lat'))
    vlon = ncfile.createVariable('lon',np.dtype('float32').char,('buoy_points','lon'))
    # results
    vu10_gefs_forecast = ncfile.createVariable('u10_gefs_forecast',np.dtype('float32').char,('buoy_points','time','ensemble_member','lat','lon'))
    vhs_gefs_forecast = ncfile.createVariable('hs_gefs_forecast',np.dtype('float32').char,('buoy_points','time','ensemble_member','lat','lon'))
    vtp_gefs_forecast = ncfile.createVariable('tp_gefs_forecast',np.dtype('float32').char,('buoy_points','time','ensemble_member','lat','lon'))
    vu10_gefs_hindcast = ncfile.createVariable('u10_gefs_hindcast',np.dtype('float32').char,('buoy_points','time','ensemble_member','lat','lon'))
    vhs_gefs_hindcast = ncfile.createVariable('hs_gefs_hindcast',np.dtype('float32').char,('buoy_points','time','ensemble_member','lat','lon'))
    vtp_gefs_hindcast = ncfile.createVariable('tp_gefs_hindcast',np.dtype('float32').char,('buoy_points','time','ensemble_member','lat','lon'))
    vu10_ndbc = ncfile.createVariable('u10_ndbc',np.dtype('float32').char,('buoy_points','time'))
    vhs_ndbc = ncfile.createVariable('hs_ndbc',np.dtype('float32').char,('buoy_points','time'))
    vtp_ndbc = ncfile.createVariable('tp_ndbc',np.dtype('float32').char,('buoy_points','time'))
    # Assign units
    vt.units = 'seconds since 1970-01-01T00:00:00+00:00'
    vblat.units = 'degrees_north' ; vblon.units = 'degrees_east'
    vlat.units = 'degrees_north' ; vlon.units = 'degrees_east'
    vu10_gefs_forecast.units = 'm/s'; vu10_gefs_hindcast.units = 'm/s'; vu10_ndbc.units = 'm/s'
    vhs_gefs_forecast.units = 'm'; vhs_gefs_hindcast.units = 'm'; vhs_ndbc.units = 'm'
    vtp_gefs_forecast.units = 's'; vtp_gefs_hindcast.units = 's'; vtp_ndbc.units = 's'
    # Allocate Data
    vt[:]=fgftime[:]; vbid[:]=np.array(ds_buoy['buoyID']).astype('str')[:]
    vblat[:]=np.array(ds_buoy['Lat']).astype('float')[:]; vblon[:]=fbmlon[:]
    vensm[:]=np.array(np.arange(0,nenm).astype('int'))[:]
    vlat[:,:]=fgflat[:,:]; vlon[:,:]=fgflon[:,:]
    vu10_gefs_forecast[:,:,:,:,:]=fgefs_forecast[:,0,:,:,:,:]; vhs_gefs_forecast[:,:,:,:,:]=fgefs_forecast[:,1,:,:,:,:]; vtp_gefs_forecast[:,:,:,:,:]=fgefs_forecast[:,2,:,:,:,:]
    vu10_gefs_hindcast[:,:,:,:,:]=fgefs_hindcast[:,0,:,:,:,:]; vhs_gefs_hindcast[:,:,:,:,:]=fgefs_hindcast[:,1,:,:,:,:]; vtp_gefs_hindcast[:,:,:,:,:]=fgefs_hindcast[:,2,:,:,:,:]
    vu10_ndbc[:,:]=fbuoy_data[:,0,:]; vhs_ndbc[:,:]=fbuoy_data[:,1,:]; vtp_ndbc[:,:]=fbuoy_data[:,2,:]
    ncfile.close()
    print(' ')
    print("Done. Netcdf ok. New file saved: "+outpath+"GEFS.BUOY.PointExtract."+date+".nc")

    stop = timeit.default_timer()
    print('Concluded in '+repr(int(round(stop - start,0)))+' seconds')
    print(" ------------------ "); print(" ")

