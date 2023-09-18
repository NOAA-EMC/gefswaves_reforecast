#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
probmaps_gefs.py

VERSION AND LAST UPDATE:
 v1.0  06/09/2023

PURPOSE:
 This program makes the probability maps based on the NOAA Global Ensemble
  Forecast System (GEFS), marine forecast. It requires the grib2 files to
  have been previously downloaded.
 The global hazards outlook (probability map) shows the probability of 
  certain variable (wind speed, significant wave height, and peak period)
  to have at least one value above certain pre-defined level within the given
  time range.
 For example, selecting the week2 time range (from day 7 to day 14), it
  calculates the probabilities of encontering at least one instant above
  qlev (defined in probmaps_gefs.yaml) during this week (for each grid point).

USAGE:
 The information is passed to the script through 4 input arguments and 
  one configuration file (probmaps_gefs.yaml).
 There are 4 input arguments:
  1) forecast cycle (YYYYMMDDHH) that will be used to compute the probabilities;
  2) initial day to define the time interval to calculate the statistics;
  3) final day to define the time interval to calculate the statistics;
  4) variable (select only one) to be processed: U10, Hs, Tp
 The configuration file saves fixed information. The most IMPORTANT variable
  to edit is the outpath, containing the location where the results (.png 
  figures) will be saved.
 This script must be run for each variable (U10, Hs, Tp) separately.
 See the probmaps_gefs.yaml for more specific information and how to calibrate
  the probabilities and customize the plots.

 Example (from linux terminal command line):
  python3 probmaps_gefs.py 2023060900 7 14 Hs
  nohup python3 probmaps_gefs.py 2023060900 7 14 Hs >> nohup_probmaps_gefs.out 2>&1 &

OUTPUT:
 .png figures (probability maps) saved in the outpath directory informed in
 the configuration file probmaps_gefs.yaml

DEPENDENCIES:
 See the imports below.

AUTHOR and DATE:
 06/09/2023: Ricardo M. Campos, first version.

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

# Pay attention to the pre-requisites and libraries
import matplotlib
matplotlib.use('Agg')
import pygrib
import xarray as xr
import matplotlib.pyplot as plt
import yaml
from matplotlib.colors import ListedColormap
from scipy.ndimage.filters import gaussian_filter
import numpy as np
import pandas as pd
import cartopy
import cartopy.crs as ccrs
# Palette and colors for plotting the figures
import matplotlib.colors as colors
import sys
import warnings; warnings.filterwarnings("ignore")
# --------------------------------------------------------------------------
sl=13 # plot style configuration
matplotlib.rcParams.update({'font.size': sl}); plt.rc('font', size=sl) 
matplotlib.rc('xtick', labelsize=sl); matplotlib.rc('ytick', labelsize=sl); matplotlib.rcParams.update({'font.size': sl})

# Input Arguments -----
# Forecast Cycle
fcycle=str(sys.argv[1])
# Forecast Lead Time (Day) and intervall
ltime1=int(sys.argv[2])
ltime2=int(sys.argv[3])
# forecast variable: U10, Hs, Tp
fvarname=str(sys.argv[4])

# Fixed configuration variables, read yaml file -----------
print(" "); print(" Reading yaml configuration file ...")
with open('probmaps_gefs.yaml', 'r') as file:
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
# Spatial percentile (grid points)
spctl=wconfig['spctl']
# gaussian filter spatial plot smoothing
gft=wconfig['gft']
# percentiles for the initial plots
pctls=np.array(wconfig['pctls']).astype('int')
# Probability levels
plevels = np.array(wconfig['plevels']).astype('float')
# Colors for the probability maps
pcolors = np.array(wconfig['pcolors']).astype('str')

# output path
if "outpath" in wconfig:
    outpath=str(wconfig['outpath'])
    if outpath[-1] != '/':
        outpath=outpath+"/"
else:
    outpath=""

umf=1. # unit conversion, when necessary

# maximum value allowed (quick quality control), variable name (grib2), and levels for the probability plot
if fvarname.upper() == "U10":
    qqvmax=wconfig['qqvmax_wnd']
    mvar=wconfig['mvar_wnd']
    qlev=np.array(wconfig['qlev_wnd']).astype('float')
    vtickd=int(wconfig['vtickd_wnd'])
    funits=str('knots')
    umf=1.94 # m/s to knots
elif fvarname.upper() == "WND":
    qqvmax=wconfig['qqvmax_wnd']
    mvar=wconfig['mvar_wnd']
    funits=str('knots')
    qlev=np.array(wconfig['qlev_wnd']).astype('float')
    vtickd=int(wconfig['vtickd_wnd'])
    umf=1.94 # m/s to knots
elif fvarname.upper() == "HS":
    qqvmax=wconfig['qqvmax_hs']
    mvar=wconfig['mvar_hs']
    funits=str('m')
    qlev=np.array(wconfig['qlev_hs']).astype('float')
    vtickd=int(wconfig['vtickd_hs'])
elif fvarname.upper() == "TP":
    qqvmax=wconfig['qqvmax_tp']
    mvar=wconfig['mvar_tp']
    funits=str('s')
    qlev=np.array(wconfig['qlev_tp']).astype('float')
    vtickd=int(wconfig['vtickd_tp'])
else:
    sys.exit(" Input variable "+fvarname+" not included in the list. Please select only one: u10, hs, tp.")

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

auxltime = np.arange((ltime1-1)*24,((ltime2)*24)+1,tres)
auxltime[auxltime>384]=384; auxltime[auxltime<0]=0

# READ WW3 Ensemble Forecast files. Appending forecast days (time intervall).
c=0
for t in range(0,auxltime.shape[0]):
    for enm in range(0,nenm):
        # fname="gefsWave."+fcycle+"/gefs.wave."+fcycle+"."+str(enm).zfill(2)+".global.0p25.f"+str(auxltime[t]).zfill(3)+".grib2"
        fname="gefsWave."+fcycle+"/gefswave"+str(enm).zfill(2)+".t"+fcycle[-2::]+"z.pgrib2f"+str(auxltime[t]).zfill(3)+".grib2"
        if c==0:
            ds = xr.open_dataset(fname, engine='cfgrib')
            wtime = np.atleast_1d(np.array(ds.time.values))
            lat = np.array(ds.latitude.values); lat = np.sort(lat); lon = np.array(ds.longitude.values)
            fmod=np.zeros((auxltime.shape[0],nenm,lat.shape[0],lon.shape[0]),'f')*np.nan
            ds.close(); del ds

        # opening and reading with pygrib (faster than xarray cfgrib)
        grbs = pygrib.open(fname); glist=list(grbs.select())
        ind=np.nan
        for i in range(0,len(glist)):
            if mvar in str(glist[i]):
                ind=int(i)
        if ind>=0:
            grb = grbs.select()[ind]
            fmod[t,enm,:,:]=np.flip(grb.values,axis=0)
        else:
            sys.exit(" Enviromental variable "+mvar+" not found in the grib2 file "+fname)

        c=c+1; del ind

    print(repr(t))

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

fmodo=np.array(fmod.reshape(auxltime.shape[0]*nenm,lat.shape[0],lon.shape[0]))

print(" 2. Initial Plots ...")
wlevels=np.linspace(0,np.nanpercentile(fmod,99.99),101)

# Percentiles
for i in range(0,pctls.shape[0]):
    plt.figure(figsize=(9,5.5))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=-90))
    ax.set_extent([slonmin+spws,slonmax-spws,slatmin+spws,slatmax-spws], crs=ccrs.PlateCarree())
    # ax.set_extent([slonmin,slonmax,slatmin,slatmax], crs=ccrs.PlateCarree())
    gl = ax.gridlines(crs=ccrs.PlateCarree(),  xlocs=range(-180,180, 20), draw_labels=True, linewidth=0.5, color='grey', alpha=0.5, linestyle='--')
    gl.xlabel_style = {'size': 9, 'color': 'k','rotation':0}; gl.ylabel_style = {'size': 9, 'color': 'k','rotation':0}
    cs=ax.contourf(lon,lat,np.nanpercentile(fmodo,pctls[i],axis=0),levels=wlevels,alpha=0.7,cmap='jet',zorder=1,extend="max",transform = ccrs.PlateCarree())
    ax.add_feature(cartopy.feature.OCEAN,facecolor=("white"))
    ax.add_feature(cartopy.feature.LAND,facecolor=("lightgrey"), edgecolor='grey',linewidth=0.5, zorder=2)
    ax.add_feature(cartopy.feature.BORDERS, edgecolor='grey', linestyle='-',linewidth=0.5, alpha=1, zorder=3)
    ax.coastlines(resolution='50m', color='dimgrey',linewidth=0.5, linestyle='-', alpha=1, zorder=4)
    title = "Percentile"+str(pctls[i]).zfill(2)+" "+fvarname+" ("+funits+"), Cycle "+fcycle[0:8]+" "+fcycle[8:10]+"Z \n"
    title += r"$\bf{"+trfi+"Valid: "+pd.to_datetime(wtime[0]+np.timedelta64(ltime1,'D')).strftime('%B %d, %Y')+" - "
    title += pd.to_datetime(wtime[0]+np.timedelta64(ltime2,'D')).strftime('%B %d, %Y')+"}$"
    ax.set_title(title); del title
    plt.tight_layout()
    ax = plt.gca(); pos = ax.get_position(); l, b, w, h = pos.bounds; cax = plt.axes([l+0.07, b-0.07, w-0.12, 0.03]) # setup colorbar axes.
    cbar = plt.colorbar(cs,cax=cax, orientation='horizontal', format='%g')
    labels = np.arange(0, wlevels.max(),vtickd).astype('int'); ticks = np.arange(0, wlevels.max(),vtickd).astype('int')
    cbar.set_ticks(ticks); cbar.set_ticklabels(labels)
    plt.axes(ax); plt.tight_layout()
    plt.savefig(outpath+"Pctl"+str(pctls[i]).zfill(2)+"_"+fvarname+"_"+fcycle+"_fcst"+str(ltime1).zfill(2)+"to"+str(ltime2).zfill(2)+"_"+ftag+".png", dpi=200, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)

    plt.close('all'); del ax
    print(" 2. Initial Plots. Percentile "+str(pctls[i]).zfill(2)+" ok.")

print(" "); print(" 3. Space-Time Cells and Probabilities ...")
# n-max expansion and reshape
fmod=np.sort(fmod,axis=0)
fmod=np.copy(fmod[-nmax::,:,:,:])
fmod=np.array(fmod.reshape(nmax*nenm,lat.shape[0],lon.shape[0]))

gspws=int(np.floor(spws/np.diff(lat).mean())/2)

probecdf=np.zeros((qlev.shape[0],lat.shape[0],lon.shape[0]),'f')
for i in range(0,qlev.shape[0]):
    for j in range(0,lat.shape[0]):
        for k in range(0,lon.shape[0]):
            if np.any(fmod[:,j,k]>0.0):
                if (j>=gspws) and (j<=lat.shape[0]-gspws) and (k>=gspws) and (k<=lon.shape[0]-gspws):
                    aux=np.array(fmod[:,(j-gspws):(j+gspws+1),:][:,:,(k-gspws):(k+gspws+1)])
                    aux=aux.reshape(nmax*nenm,aux.shape[1]*aux.shape[2])
                    aux=np.sort(aux,axis=1)
                    ind=np.where(np.mean(aux,axis=0)>=0.)
                    if np.size(ind)>0:
                        aux=np.array(aux[:,ind[0]])            
                        aux=np.array(aux[:,int(np.floor(aux.shape[1]*(spctl/100)))::])
                        aux=aux.reshape(aux.shape[0]*aux.shape[1])
                        # 1054 for spws=2.0 and nmax=2 and spctl=80 (one week)
                        # 682 for spws=2.0 and nmax=2 and spctl=87 (one week)
                        # 558 for spws=2.0 and nmax=2 and spctl=90 (one week)
                        probecdf[i,j,k] = np.size(aux[aux>qlev[i]]) / np.size(aux[aux>0.])

                    del aux

print(" 3. Space-Time Cells and Probabilities ... OK"); print(" ")

# PLOTS
print(" 4. Probability Maps ...")
# Probability levels and labels for the plots
clabels=[];clevels=[]
for i in range(0,np.size(plevels)-1):
    clabels=np.append(clabels,">"+str(int(plevels[i]*100)).zfill(2)+"%")
    clevels=np.append(clevels,(plevels[i]+plevels[i+1])/2)

cmap = ListedColormap(pcolors)

for i in range(0,qlev.shape[0]):
    plt.figure(figsize=(9,5.5))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=-90))
    ax.set_extent([slonmin+spws,slonmax-spws,slatmin+spws,slatmax-spws], crs=ccrs.PlateCarree())
    # ax.set_extent([slonmin,slonmax,slatmin,slatmax], crs=ccrs.PlateCarree())
    gl = ax.gridlines(crs=ccrs.PlateCarree(),  xlocs=range(-180,180, 20), draw_labels=True, linewidth=0.5, color='grey', alpha=0.5, linestyle='--')
    gl.xlabel_style = {'size': 9, 'color': 'k','rotation':0}; gl.ylabel_style = {'size': 9, 'color': 'k','rotation':0}
    cs=ax.contourf(lon,lat,gaussian_filter(probecdf[i,:,:],gft),levels=plevels,alpha=0.7,cmap=cmap,zorder=1,transform = ccrs.PlateCarree())
    ax.add_feature(cartopy.feature.OCEAN,facecolor=("white"))
    ax.add_feature(cartopy.feature.LAND,facecolor=("lightgrey"), edgecolor='grey',linewidth=0.5, zorder=2)
    ax.add_feature(cartopy.feature.BORDERS, edgecolor='grey', linestyle='-',linewidth=0.5, alpha=1, zorder=3)
    ax.coastlines(resolution='50m', color='dimgrey',linewidth=0.5, linestyle='-', alpha=1, zorder=4)
    title = "Prob "+fvarname+">"+str(qlev[i]).zfill(1)+funits+", Cycle "+fcycle[0:8]+" "+fcycle[8:10]+"Z \n"
    title += r"$\bf{"+trfi+"Valid: "+pd.to_datetime(wtime[0]+np.timedelta64(ltime1,'D')).strftime('%B %d, %Y')+" - "
    title += pd.to_datetime(wtime[0]+np.timedelta64(ltime2,'D')).strftime('%B %d, %Y')+"}$"
    ax.set_title(title); del title
    plt.tight_layout()
    ax = plt.gca(); pos = ax.get_position(); l, b, w, h = pos.bounds; cax = plt.axes([l+0.07, b-0.07, w-0.12, 0.03]) # setup colorbar axes.
    cbar = plt.colorbar(cs,cax=cax, orientation='horizontal',ticks=clevels, format='%g')
    cbar.ax.set_xticklabels(clabels)
    cbar.ax.tick_params(length=0)
    for label in cbar.ax.get_xticklabels():
        label.set_weight('bold')

    plt.axes(ax); plt.tight_layout()
    plt.savefig(outpath+"ProbMap_"+fvarname+"_"+str(qlev[i]).zfill(1)+"_"+fcycle+"_fcst"+str(ltime1).zfill(2)+"to"+str(ltime2).zfill(2)+"_"+ftag+".png", dpi=200, facecolor='w', edgecolor='w',
            orientation='portrait', papertype=None, format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)

    plt.close('all'); del ax
    print("   Plot ... qlev "+repr(qlev[i]))

print(" 4. Probability Plots ... OK"); print(" ")

