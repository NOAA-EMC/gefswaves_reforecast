#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
hindcast_reference_archives.py

VERSION AND LAST UPDATE:
 v1.0  06/09/2023
 v1.1  05/09/2024

PURPOSE:
 It creates and plots the assembled hindcast from GEFSv12 forecast archives,
  to build a reference, served as the ground truth to validate the
  probability maps (generated with probmaps_gefs_refcst.py).
 The code does not read directly the .grib2 operational files but the 
  netcdf files post-processed by postproc_AWSarchive_GEFS_Field.sh and
  jobmanager_postproc_AWSarchive_GEFS_Field.sh. 
 The plot colors, pallete, and contours were conceived to highlight the
  intensity levels used in probmaps_gefs_refcst.py and defined in 
  probmaps_gefs_refcst.py

USAGE:
 This code mimics the main program probmaps_gefs_refcst.py to be consistent
  with the probability maps and allow a direct comparison. Therefore the
  inputs are the same.
 The information is passed to the script through 6 input arguments and 
  one configuration file (probmaps_gefs.yaml).
 There are 4 input arguments:
  1) path and name of .yaml descriptive file (probmaps_gefs_refcst.yaml);
  2) forecast cycle (YYYYMMDDHH) that will be used to compute the probabilities;
  3) initial day to define the time interval to calculate the statistics;
  4) final day to define the time interval to calculate the statistics;
  5) variable (select only one) to be processed: U10, Hs;
  6) name of the event, which will be included in the output name of the figure;

 Make sure the information and paths in probmaps_gefs_refcst.yaml are correct.
 This code is run inside probmaps_gefs_refcst_fromList.sh which reads a list
  of important events to be re-analyzed list_events.txt

OUTPUT:
 .png figure (HindcastReference_*.png) saved in the outpath directory 
  informed in the configuration file probmaps_gefs_refcst.yaml

DEPENDENCIES:
 See the imports below.

AUTHOR and DATE:
 06/09/2023: Ricardo M. Campos, first version named probmaps_reference_archives.py
 05/09/2024: Ricardo M. Campos, improvements in the plot and colors.

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

import matplotlib
matplotlib.use('Agg')
import xarray as xr
import matplotlib.pyplot as plt
import yaml
import numpy as np
import pandas as pd
from calendar import timegm
import time
from time import strptime
import cartopy
import cartopy.crs as ccrs
import sys
import gc
import warnings; warnings.filterwarnings("ignore")
# --------------------------------------------------------------------------
sl=13 # plot style configuration
matplotlib.rcParams.update({'font.size': sl}); plt.rc('font', size=sl) 
matplotlib.rc('xtick', labelsize=sl); matplotlib.rc('ytick', labelsize=sl); matplotlib.rcParams.update({'font.size': sl})

def build_gefs_hindcast(date,mvar,auxftime,gefspath):

    # cycle time
    ctime=np.array(np.arange(auxftime.min(),auxftime.max()+1,24)*3600.+timegm( strptime(date,'%Y%m%d') )).astype('double')
    auxltime=np.array(np.arange(0,24,6)).astype('int'); auxltimef=np.array([0]).astype('int')

    if str.lower(mvar)=='hs':
        mvar=str('HTSGW_surface')
    if str.lower(mvar)=='wsp' or str.lower(mvar)=='u10':
        mvar=str('WIND_surface')

    c=0; tt=0
    for ct in range(0,len(ctime)):

        if ct==int(len(ctime)-1):
            auxlt=auxltimef
        else:
            auxlt=auxltime

        for enm in range(0,nenm):
            cdate=str(time.gmtime(ctime[ct])[0])+str(time.gmtime(ctime[ct])[1]).zfill(2)+str(time.gmtime(ctime[ct])[2]).zfill(2)
            fname=gefspath+"gefs.wave."+cdate+"."+str(int(enm)).zfill(2)+".global.0p25.nc"
            print(" build_gefs_hindcast "+fname)

            if c==0:
                ds = xr.open_dataset(fname)
                wtime = np.atleast_1d(np.array(ds.time.values))
                lat = np.array(ds.latitude.values); lat = np.sort(lat); lon = np.array(ds.longitude.values)
                fmod=np.zeros((len(auxftime),nenm,len(lat),len(lon)),'f')*np.nan
                ds.close(); del ds

            ds = xr.open_dataset(fname)
            fmod[tt:tt+len(auxlt),enm,:,:]=np.array(ds[mvar].values[0:len(auxlt),:,:])
            ds.close()
            c=c+1

        tt=tt+len(auxlt)

    # wtime = np.array( auxftime*3600 + (wtime.astype('int64') // 10**9) ).astype('double')
    return fmod, lat, lon, wtime

if __name__ == "__main__":

    # Input Arguments -----
    # .yaml configuration file
    fconfig=str(sys.argv[1])
    # Forecast Cycle
    fcycle=str(sys.argv[2]); fcdate=str(fcycle[0:8])
    # Forecast Lead Time (Day) and intervall
    ltime1=int(sys.argv[3])
    ltime2=int(sys.argv[4])
    # forecast variable: WS10, Hs
    fvarname=str(sys.argv[5])
    if len(sys.argv) >= 5:
        name_event="_"+str(sys.argv[6])
    else:
        name_event=str('')

    # Fixed configuration variables, read yaml file -----------
    print(" "); print(" Reading yaml configuration file ...")
    with open(fconfig, 'r') as file:
        wconfig = yaml.safe_load(file)

    ftag=str(wconfig['ftag'])

    # number of ensemble members
    nenm=wconfig['nenm']
    # time resolution
    tres=wconfig['tres']
    # n-max expansion
    nmax=wconfig['nmax']
    # spatial window size (diameter), square (degrees)
    spws=wconfig['spws']
    # Plotting Area
    slonmin=wconfig['lonmin']-spws; slonmax=wconfig['lonmax']+spws
    slatmin=wconfig['latmin']-spws; slatmax=wconfig['latmax']+spws

    # Path of GEFSv12 files
    gefspath=wconfig['gefspath']
    if gefspath[-1] != '/':
        gefspath=gefspath+"/"

    # output path
    outpath=str(wconfig['outpath'])
    if outpath[-1] != '/':
        outpath=outpath+"/"

    umf=1. # unit conversion, when necessary

    # maximum value allowed (quick quality control), variable name (grib2), and levels for the probability plot
    if fvarname.upper() == "WS10" or fvarname.upper() == "WND" or fvarname.upper() == "U10":
        qqvmax=wconfig['qqvmax_wnd']
        mvar=wconfig['mvar_wnd']
        qlev=np.array(wconfig['qlev_wnd']).astype('float')
        vtickd=int(wconfig['vtickd_wnd'])
        funits=str('knots')
        umf=1.94 # m/s to knots
    elif fvarname.upper() == "HS":
        qqvmax=wconfig['qqvmax_hs']
        mvar=wconfig['mvar_hs']
        funits=str('m')
        qlev=np.array(wconfig['qlev_hs']).astype('float')
        vtickd=int(wconfig['vtickd_hs'])
    else:
        sys.exit(" Input variable "+fvarname+" not included in the list. Please select only one: WS10, Hs.")

    del wconfig
    print(" Reading yaml configuration file, OK."); print(" ")

    # Time range forecast intervall string, for the plots
    if (ltime2-ltime1)<0:
        aux=np.copy(ltime1); ltime1=np.copy(ltime2)
        ltime2=np.copy(aux); del aux

    if ltime1==1 and (ltime2-ltime1)==7:
        trfi=str("Week 1 - ")
    elif ltime1==7 and (ltime2-ltime1)==7:
        trfi=str("Week 2 - ")
    elif ltime1==14 and (ltime2-ltime1)==7:
        trfi=str("Week 3 - ")
    elif ltime1==21 and (ltime2-ltime1)==7:
        trfi=str("Week 4 - ")
    elif ltime1==28 and (ltime2-ltime1)==7:
        trfi=str("Week 5 - ")
    elif (ltime2-ltime1)>0:
        trfi=str("Days "+str(ltime1)+"-"+str(ltime2)+" , ")
    elif ltime2==ltime1:
        trfi=str("Day "+str(ltime1)+" , ")
    else:
        trfi=''

    print(" "); print(" 1. Reading Forecast Data ...")

    # Read GEFS and build 24-h slices (Hindcast)
    auxftime = np.arange((ltime1-1)*24,((ltime2)*24)+1,tres)
    auxftime[auxftime>384]=384; auxftime[auxftime<0]=0

    fmod, lat, lon, wtime = build_gefs_hindcast(fcdate,mvar,auxftime,gefspath)
    # -----

    # Quick simple quality control
    fmod[fmod>=qqvmax]=np.nan; fmod[fmod<0.]=np.nan
    # Unit conversion
    fmod=fmod*umf
    # Select domain of interest
    indlat=np.where((lat>=(slatmin-spws))&(lat<=(slatmax+spws)))
    indlon=np.where((lon>=(slonmin-spws))&(lon<=(slonmax+spws)))
    if np.size(indlon)>0 and np.size(indlat)>0:
        lat=np.copy(lat[indlat[0]]); lon=np.copy(lon[indlon[0]])
        fmod=np.copy(fmod[:,:,indlat[0],:][:,:,:,indlon[0]])
    else:
        sys.exit(" Min/Max lat and lon incorrect when applied to the forecast file. Check longitude standards.")

    del indlat, indlon
    print(" 1. Forecast Data ... OK"); print(" ")

    print(" 2. Plots ...")

    wdata = np.nanmean(np.sort(np.nanmean(fmod,axis=1),axis=0)[-nmax::,:,:],axis=0)

    blevels=np.array(qlev)
    if fvarname=="Hs":
        plevels=np.array(np.arange(blevels[0],blevels[-1]+1,1.0)).astype('int')
    else:
        plevels=np.array(np.arange(blevels[0],blevels[-1]+1,2.0)).astype('int')

    wlevels=np.arange(blevels[0],blevels[-1]*1.01,0.1)

    plt.figure(figsize=(9,5.5))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=-90))
    ax.set_extent([slonmin+spws,slonmax-spws,slatmin+spws,slatmax-spws], crs=ccrs.PlateCarree())
    # ax.set_extent([slonmin,slonmax,slatmin,slatmax], crs=ccrs.PlateCarree())
    gl = ax.gridlines(crs=ccrs.PlateCarree(),  xlocs=range(-180,180, 20), draw_labels=True, linewidth=0.5, color='grey', alpha=0.5, linestyle='--')
    gl.xlabel_style = {'size': 9, 'color': 'k','rotation':0}; gl.ylabel_style = {'size': 9, 'color': 'k','rotation':0}
    cs=ax.contourf(lon,lat,wdata,levels=wlevels,vmin=blevels[0],alpha=0.7,cmap='jet',zorder=1,extend="max",transform = ccrs.PlateCarree())
    ct=ax.contour(lon,lat,wdata,levels=blevels,colors='dimgrey', linewidths=1.,transform = ccrs.PlateCarree())
    plt.clabel(ct, inline=True, fontsize=6,colors='black',fmt = '%.0f')
    ax.add_feature(cartopy.feature.OCEAN,facecolor=("white"))
    ax.add_feature(cartopy.feature.LAND,facecolor=("lightgrey"), edgecolor='grey',linewidth=0.5, zorder=2)
    ax.add_feature(cartopy.feature.BORDERS, edgecolor='grey', linestyle='-',linewidth=0.5, alpha=1, zorder=3)
    ax.coastlines(resolution='50m', color='dimgrey',linewidth=0.5, linestyle='-', alpha=1, zorder=4)
    title = "Reference  "+fvarname+" ("+funits+") \n"
    title += r"$\bf{Valid: "+pd.to_datetime(wtime[0]+np.timedelta64(1,'D')).strftime('%B %d, %Y')+" - "
    title += pd.to_datetime(wtime[0]+np.timedelta64(1+ltime2-ltime1,'D')).strftime('%B %d, %Y')+"}$"
    ax.set_title(title); del title
    plt.tight_layout()
    ax = plt.gca(); pos = ax.get_position(); l, b, w, h = pos.bounds; cax = plt.axes([l+0.07, b-0.07, w-0.12, 0.03]) # setup colorbar axes.
    cbar = plt.colorbar(cs,cax=cax, orientation='horizontal', format='%g')
    labels = np.array(plevels).astype('int'); ticks = np.array(plevels).astype('int')
    cbar.set_ticks(ticks); cbar.set_ticklabels(labels)
    plt.axes(ax); plt.tight_layout()

    figname = outpath+"HindcastReference_"+fvarname+"_"+fcdate+"_fcst"+str(ltime1).zfill(2)+"to"+str(ltime2).zfill(2)+name_event
    plt.savefig(figname+".png", dpi=200, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)

    plt.close('all'); del ax
    print(" Plots Ok.")
    gc.collect()

