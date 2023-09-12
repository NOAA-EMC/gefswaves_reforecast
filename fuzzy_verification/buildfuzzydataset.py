#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
buildfuzzydataset.py

VERSION AND LAST UPDATE:
 v1.0  09/10/2023

PURPOSE:
 Python script to build a dataset with GEFS forecast, GEFS forecast (hindcasted 24-h slices),
  GDAS, and NDBC buoy data. It uses the NDBC buoy point positions to select the nearest model point and
  neighbouring points using a 10X10 degree window (wl=5) centered at the NDBC buoy point.
  This is intended for fuzzy verification.

USAGE:
 This code is run for one forecast cycle only. A shell script can execute this script multiple times.
 python3 buildfuzzydataset.py 20220701

OUTPUT:
 Netcdf file GEFS.GDAS.BUOY.3PointExtract.YYYYMMDD.nc saved in the given outpath.

DEPENDENCIES:
 See setup.py and the imports below.

AUTHOR and DATE:
 09/10/2023: Ricardo M. Campos, first version.

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

import matplotlib
matplotlib.use('Agg')
from matplotlib.mlab import *
import pygrib
import numpy as np
import sys
import pandas as pd
import xarray as xr
import netCDF4 as nc
import time
import timeit
from time import strptime
from calendar import timegm
from datetime import datetime, timedelta
import wproc
import wread
import warnings; warnings.filterwarnings("ignore")
# netcdf format
fnetcdf="NETCDF4"


def read_gefs(date,blat,blon,wl,path_gefs):
    # variable names read by pygrib
    mvar=np.array(["Wind speed","Significant height of combined wind waves and swell","Primary wave mean period"]).astype('str')

    # buoy lon standard
    if blon<0:
        blon=blon+360

    c=0
    auxltime=np.array(np.arange(0,384+1,6)).astype('int')
    for t in range(0,auxltime.shape[0]):
        for enm in range(0,nenm):

            fname=path_gefs+"GEFSv12Waves_"+date+"/gefs.wave."+date+"."+str(int(enm)).zfill(2)+".global.0p25.f"+str(int(auxltime[t])).zfill(3)+".grib2"
            print(" read_gefs "+fname)
            if c==0:
                ds = xr.open_dataset(fname, engine='cfgrib')
                wtime = np.atleast_1d(np.array(ds.time.values))
                lat = np.array(ds.latitude.values); lat = np.sort(lat); lon = np.array(ds.longitude.values)

                dlat=int(np.round(wl/np.mean(np.diff(lat))))
                indlat=np.min(np.where(np.abs(lat-blat) == np.min(np.abs(lat-blat))))
                indlat=np.array(np.arange(int(indlat-dlat),int(indlat+dlat+1))).astype('int')

                dlon=int(np.round(wl/np.mean(np.diff(lon))))
                indlon=np.min(np.where(np.abs(lon-blon) == np.min(np.abs(lon-blon))))
                indlon=np.array(np.arange(int(indlon-dlon),int(indlon+dlon+1))).astype('int')

                fmod=np.zeros((len(mvar),len(auxltime),nenm,len(indlat),len(indlon)),'f')*np.nan
                ds.close(); del ds

            # opening and reading with pygrib (faster than xarray cfgrib)
            grbs = pygrib.open(fname); glist=list(grbs.select())
            for i in range(0,len(mvar)):
                ind=np.nan
                for j in range(0,len(glist)):
                    if mvar[i] in str(glist[j]):
                        ind=int(j)

                if ind>=0:
                    grb = grbs.select()[ind]
                    fmod[i,t,enm,:,:]=np.array(np.flip(grb.values,axis=0))[indlat,indlon]

                del ind

            grbs.close(); del grbs
            c=c+1
            # print(repr(c)); print(repr(t)+" "+repr(enm))

    wtime = np.array( auxltime*3600 + (wtime.astype('int64') // 10**9) ).astype('double')
    return fmod, np.array(lat[indlat]), np.array(lon[indlon]), wtime


def build_gefs_hindcast(date,blat,blon,wl,path_gefs):

    # variable names read by pygrib
    mvar=np.array(["Wind speed","Significant height of combined wind waves and swell","Primary wave mean period"]).astype('str')
    # cycle time
    ctime=np.array(np.arange(0,384+1,24)*3600.+timegm( strptime(date,'%Y%m%d') )).astype('double')

    if blon<0:
        blon=blon+360

    auxftime=np.array(np.arange(0,384+1,6)).astype('int')
    auxltime=np.array(np.arange(0,24,6)).astype('int'); auxltimef=np.array([0]).astype('int')

    c=0; tt=0
    for ct in range(0,len(ctime)):

        if ct==16:
            auxlt=auxltimef
        else:
            auxlt=auxltime

        for t in range(0,auxlt.shape[0]):
            for enm in range(0,nenm):
                cdate=str(time.gmtime(ctime[ct])[0])+str(time.gmtime(ctime[ct])[1]).zfill(2)+str(time.gmtime(ctime[ct])[2]).zfill(2)
                fname=path_gefs+"GEFSv12Waves_"+cdate+"/gefs.wave."+cdate+"."+str(int(enm)).zfill(2)+".global.0p25.f"+str(int(auxlt[t])).zfill(3)+".grib2"
                print(" build_gefs_hindcast "+fname)
                if c==0:
                    ds = xr.open_dataset(fname, engine='cfgrib')
                    wtime = np.atleast_1d(np.array(ds.time.values))
                    lat = np.array(ds.latitude.values); lat = np.sort(lat); lon = np.array(ds.longitude.values)

                    dlat=int(np.round(wl/np.mean(np.diff(lat))))
                    indlat=np.min(np.where(np.abs(lat-blat) == np.min(np.abs(lat-blat))))
                    indlat=np.array(np.arange(int(indlat-dlat),int(indlat+dlat+1))).astype('int')

                    dlon=int(np.round(wl/np.mean(np.diff(lon))))
                    indlon=np.min(np.where(np.abs(lon-blon) == np.min(np.abs(lon-blon))))
                    indlon=np.array(np.arange(int(indlon-dlon),int(indlon+dlon+1))).astype('int')

                    fmod=np.zeros((len(mvar),len(auxftime),nenm,len(indlat),len(indlon)),'f')*np.nan
                    ds.close(); del ds

                # opening and reading with pygrib (faster than xarray cfgrib)
                grbs = pygrib.open(fname); glist=list(grbs.select())
                for i in range(0,len(mvar)):
                    ind=np.nan
                    for j in range(0,len(glist)):
                        if mvar[i] in str(glist[j]):
                            ind=int(j)

                    if ind>=0:
                        grb = grbs.select()[ind]
                        fmod[i,tt,enm,:,:]=np.array(np.flip(grb.values,axis=0))[indlat,indlon]

                    del ind

                grbs.close(); del grbs
                c=c+1
                # print(repr(c)); print(repr(t)+" "+repr(enm)); print(repr(tt)); print(" ")

            tt=tt+1

    wtime = np.array( auxftime*3600 + (wtime.astype('int64') // 10**9) ).astype('double')
    return fmod, np.array(lat[indlat]), np.array(lon[indlon]), wtime


def read_gdas(date,blat,blon,wl,path_gdas):

    # variable names
    mvar=np.array(['WIND_surface','HTSGW_surface','PERPW_surface']).astype('str')

    # buoy lon standard
    if blon<0:
        blon=blon+360

    auxltime=np.array(np.arange(0,384+1,6)).astype('int')
    atime=np.array(auxltime*3600.+timegm( strptime(date,'%Y%m%d') )).astype('double')
    # month and year
    datestrg1=str(time.gmtime(atime[0])[0])+str(time.gmtime(atime[0])[1]).zfill(2)
    datestrg2=str(time.gmtime(atime[-1])[0])+str(time.gmtime(atime[-1])[1]).zfill(2)

    print(" read_gdas "+str(date))
    if datestrg1 == datestrg2 :
        f=nc.Dataset(path_gdas+"gdaswave."+str(datestrg1)+".global.0p25.nc")
        t=f.variables['time'][:]; indt=np.where(np.isin(t,atime)==True)[0]

        lat=f.variables['latitude'][:]; lon=f.variables['longitude'][:]

        dlat=int(np.round(wl/np.mean(np.diff(lat))))
        indlat=np.min(np.where(np.abs(lat-blat) == np.min(np.abs(lat-blat))))
        indlat=np.array(np.arange(int(indlat-dlat),int(indlat+dlat+1))).astype('int')

        dlon=int(np.round(wl/np.mean(np.diff(lon))))
        indlon=np.min(np.where(np.abs(lon-blon) == np.min(np.abs(lon-blon))))
        indlon=np.array(np.arange(int(indlon-dlon),int(indlon+dlon+1))).astype('int')  

        fmod=np.zeros((len(mvar),len(atime),len(indlat),len(indlon)),'f')*np.nan

        for i in range(0,len(mvar)):
            fmod[i,:,:,:]=np.array(f.variables[mvar[i]][indt,indlat,indlon])

    else:

        f1=nc.Dataset(path_gdas+"gdaswave."+str(datestrg1)+".global.0p25.nc")
        f2=nc.Dataset(path_gdas+"gdaswave."+str(datestrg2)+".global.0p25.nc")

        t=np.append(f1.variables['time'][:],f2.variables['time'][:]) ; indt=np.where(np.isin(t,atime)==True)[0]

        lat=f1.variables['latitude'][:]; lon=f1.variables['longitude'][:]

        dlat=int(np.round(wl/np.mean(np.diff(lat))))
        indlat=np.min(np.where(np.abs(lat-blat) == np.min(np.abs(lat-blat))))
        indlat=np.array(np.arange(int(indlat-dlat),int(indlat+dlat+1))).astype('int')

        dlon=int(np.round(wl/np.mean(np.diff(lon))))
        indlon=np.min(np.where(np.abs(lon-blon) == np.min(np.abs(lon-blon))))
        indlon=np.array(np.arange(int(indlon-dlon),int(indlon+dlon+1))).astype('int') 

        fmod=np.zeros((len(mvar),len(atime),len(indlat),len(indlon)),'f')*np.nan

        for i in range(0,len(mvar)):
            fmod[i,:,:,:]=np.array(np.append(f1.variables[mvar[i]][:,indlat,indlon],f1.variables[mvar[i]][:,indlat,indlon],axis=0)[indt,:,:])

    return fmod, np.array(lat[indlat]), np.array(lon[indlon]), atime


def read_buoy(date,bid,anh,path_buoy):

    auxltime=np.array(np.arange(0,384+1,6)).astype('int')
    atime=np.array(auxltime*3600.+timegm( strptime(date,'%Y%m%d') )).astype('double')
    # year
    datestrg1=str(time.gmtime(atime[0])[0]); datestrg2=str(time.gmtime(atime[-1])[0])
    print(" read_buoy "+str(bid)+" "+str(date))
    if datestrg1 == datestrg2 :
        bdata=wread.tseriesnc_ndbc(path_buoy+str(bid)+"h"+datestrg1+".nc",anh=anh)
        bt=np.array(bdata['time']).astype('double')
        fbdata=np.array([bdata['wind_spd']])
        fbdata=np.append(fbdata,[bdata['hs']],axis=0)
        fbdata=np.append(fbdata,[bdata['tp']],axis=0)
        ndbclat=float(bdata['latitude']); ndbclon=float(bdata['longitude'])
        del bdata
    else:
        bdata1=wread.tseriesnc_ndbc(path_buoy+str(bid)+"h"+datestrg1+".nc")
        bdata2=wread.tseriesnc_ndbc(path_buoy+str(bid)+"h"+datestrg2+".nc")
        bt=np.array(np.append(bdata1['time'],bdata2['time'])).astype('double')
        fbdata=np.array([np.append(bdata1['wind_spd'],bdata2['wind_spd'])])
        fbdata=np.append(fbdata,[np.append(bdata1['hs'],bdata2['hs'])],axis=0)
        fbdata=np.append(fbdata,[np.append(bdata1['tp'],bdata2['tp'])],axis=0)
        ndbclat=float(bdata1['latitude']); ndbclon=float(bdata1['longitude'])
        del bdata1,bdata2

    fdata=np.zeros((3,len(atime)),'f')*np.nan
    for i in range(0,len(atime)):
        indt=np.where(np.abs(atime[i]-bt)<3600)
        if np.size(indt)>0:
            indt=indt[0]
            for j in range(0,fbdata.shape[0]):
                fdata[j,i]=np.nanmean(fbdata[j,indt])

            del indt

    return fdata, ndbclat, ndbclon, atime


if __name__ == "__main__":

    # Input argument (cycle time only):
    date=str(sys.argv[1])

    # start time
    start = timeit.default_timer()

    # spatial window 
    wl=5 # (°)
    # number of ensemble members
    nenm=int(30+1)
    # Data paths
    path_gefs="/work/noaa/marine/ricardo.campos/data/archiveOPruns/GEFSv12Waves_AWS/" # GEFSv12Waves_20220720 gefs.wave.20220722.20.global.0p25.f084.grib2
    path_gdas="/work/noaa/marine/ricardo.campos/data/archiveOPruns/GDASwave/" # gdaswave.202206.global.0p25.nc
    path_buoy="/work/noaa/marine/ricardo.campos/data/buoys/NDBC/ncformat/wparam/" # 41048h2022.nc
    outpath="/work/noaa/marine/ricardo.campos/work/analysis/3assessments/fuzzy_verification/data/"

    # Read Buoy points
    ds_buoy = pd.read_csv('/work/noaa/marine/ricardo.campos/work/analysis/3assessments/fuzzy_verification/buoys_sel.txt',comment='#',delimiter=r"\s+")

    for i in range(0,len(ds_buoy)):
        bid=ds_buoy['buoyID'][i]; anh=ds_buoy['AnH'][i]
        blat=float(ds_buoy['Lat'][i]); blon=float(ds_buoy['Lon'][i])

        # GEFS forecast
        gefs_forecast, gflat, gflon, gftime = read_gefs(date,blat,blon,wl,path_gefs)
        # GEFS sliced hindcasted
        gefs_hindcast, ghlat, ghlon, ghtime = build_gefs_hindcast(date,blat,blon,wl,path_gefs)
        # GDAS
        gdas, gdlat, gdlon, gdtime = read_gdas(date,blat,blon,wl,path_gdas)
        # Buoy data
        buoy_data, bdlat, bdlon, bdtime = read_buoy(date,bid,anh,path_buoy)

        if (np.array_equal(gflat,ghlat)==True) and (np.array_equal(gflat,gdlat)==True):
            if (np.array_equal(gftime,ghtime)==True) and (np.array_equal(gftime,gdtime)==True):
                if i==0:
                    fgefs_forecast = np.array([gefs_forecast])
                    fgefs_hindcast = np.array([gefs_hindcast])
                    fgdas = np.array([gdas])
                    fbuoy_data = np.array([buoy_data])
                    fgftime = np.array(gftime)
                    fgflat = np.array([gflat])
                    fgflon = np.array([gflon])
                else:
                    fgefs_forecast = np.append(fgefs_forecast,np.array([gefs_forecast]),axis=0)
                    fgefs_hindcast = np.append(fgefs_hindcast,np.array([gefs_hindcast]),axis=0)
                    fgdas = np.append(fgdas,np.array([gdas]),axis=0)
                    fbuoy_data = np.append(fbuoy_data,np.array([buoy_data]),axis=0)
                    fgflat = np.append(fgflat,np.array([gflat]),axis=0)
                    fgflon = np.append(fgflon,np.array([gflon]),axis=0)

        del gefs_forecast,gefs_hindcast,gdas,buoy_data,gftime,gflat,gflon
        print(" "); print(" --- Point "+repr(bid)); print(" ")

    # keep model lon same standard as NDBC lon (-180 to 180)
    fbmlon=np.array(ds_buoy['Lon']).astype('float'); fbmlon[fbmlon>180]=fbmlon[fbmlon>180]-360.

    # Save netcdf output file
    ncfile = nc.Dataset(outpath+"GEFS.GDAS.BUOY.3PointExtract."+date+".nc", "w", format=fnetcdf)
    ncfile.history="Data extracted for 3 NDBC buoy position and neighbouring points using 10X10 degree window."
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
    vu10_gdas = ncfile.createVariable('u10_gdas',np.dtype('float32').char,('buoy_points','time','lat','lon'))
    vhs_gdas = ncfile.createVariable('hs_gdas',np.dtype('float32').char,('buoy_points','time','lat','lon'))
    vtp_gdas = ncfile.createVariable('tp_gdas',np.dtype('float32').char,('buoy_points','time','lat','lon'))
    vu10_ndbc = ncfile.createVariable('u10_ndbc',np.dtype('float32').char,('buoy_points','time'))
    vhs_ndbc = ncfile.createVariable('hs_ndbc',np.dtype('float32').char,('buoy_points','time'))
    vtp_ndbc = ncfile.createVariable('tp_ndbc',np.dtype('float32').char,('buoy_points','time'))
    # Assign units
    vt.units = 'seconds since 1970-01-01T00:00:00+00:00'
    vblat.units = 'degrees_north' ; vblon.units = 'degrees_east'
    vlat.units = 'degrees_north' ; vlon.units = 'degrees_east'
    vu10_gefs_forecast.units = 'm/s'; vu10_gefs_hindcast.units = 'm/s'; vu10_gdas.units = 'm/s'; vu10_ndbc.units = 'm/s'
    vhs_gefs_forecast.units = 'm'; vhs_gefs_hindcast.units = 'm'; vhs_gdas.units = 'm'; vhs_ndbc.units = 'm'
    vtp_gefs_forecast.units = 's'; vtp_gefs_hindcast.units = 's'; vtp_gdas.units = 's'; vtp_ndbc.units = 's'
    # Allocate Data
    vt[:]=fgftime[:]; vbid[:]=np.array(ds_buoy['buoyID']).astype('str')[:]
    vblat[:]=np.array(ds_buoy['Lat']).astype('float')[:]; vblon[:]=fbmlon[:]
    vensm[:]=np.array(np.arange(0,nenm).astype('int'))[:]
    vlat[:,:]=fgflat[:,:]; vlon[:,:]=fgflon[:,:]
    vu10_gefs_forecast[:,:,:,:,:]=fgefs_forecast[:,0,:,:,:,:]; vhs_gefs_forecast[:,:,:,:,:]=fgefs_forecast[:,1,:,:,:,:]; vtp_gefs_forecast[:,:,:,:,:]=fgefs_forecast[:,2,:,:,:,:]
    vu10_gefs_hindcast[:,:,:,:,:]=fgefs_hindcast[:,0,:,:,:,:]; vhs_gefs_hindcast[:,:,:,:,:]=fgefs_hindcast[:,1,:,:,:,:]; vtp_gefs_hindcast[:,:,:,:,:]=fgefs_hindcast[:,2,:,:,:,:]
    vu10_gdas[:,:,:,:]=fgdas[:,0,:,:,:]; vhs_gdas[:,:,:,:]=fgdas[:,1,:,:,:]; vtp_gdas[:,:,:,:]=fgdas[:,2,:,:,:]
    vu10_ndbc[:,:]=fbuoy_data[:,0,:]; vhs_ndbc[:,:]=fbuoy_data[:,1,:]; vtp_ndbc[:,:]=fbuoy_data[:,2,:]
    ncfile.close()
    print(' ')
    print("Done. Netcdf ok. New file saved: "+outpath+"GEFS.GDAS.BUOY.3PointExtract."+date+".nc")

    stop = timeit.default_timer()
    print('Concluded in '+repr(int(round(stop - start,0)))+' seconds')


