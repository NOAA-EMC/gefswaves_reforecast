#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
checkww3result.py

VERSION AND LAST UPDATE:
 v1.0  05/16/2022

PURPOSE:
 Build a panel with multiple plots to check WAVEWATCHIII global 
  simulation was run correctly. It includes field outputs, table with
  time-series and spectral output.
 The verification shows the WAVEWATCHIII results and also includes a
  comparison with a NDBC buoy to ensure everything is correct.

USAGE:
 Three WAVEWATCHIII output types must exist: ww3gefs.*_field.grib2,
  ww3gefs.*_spec.nc, and ww3gefs.*_tab.nc, which will be read by this
  program. Additionally, NDBC buoys must be available (netcdf format),
  obtained with retrieve_ndbc_nc.py.
 Two input arguments must be informed:
  - date (YYYYMMDD)
  - tag, which is a string that will be included in the output file name
 Moreover, two paths must be informed: wrun and pbuoys (see below).
 See also https://github.com/NOAA-EMC/WW3-tools

OUTPUT:
 One png figure with multiple subplots: wave fiels of wind speed, Hs, and
  Tm; directional spectrum; and time-series of wind speed, Hs, and Tm, 
  where WAVEWATCHIII is compared with NDBC buoys.
 The panel/figure includes four columns and three lines: first line is 
  the first time step (to verify ww3 read the restart properly and the
  wind is ok), and the second line is the last time step (to confirm it
  went well through the entire run). The third line shows the time series,
  comparing with the NDBC buoy available (it loops through the data 
  searching for a buoy with valid numbers). Finally, the information on 
  the right corner is a general view to confirm it is ok. 
 The goal is to run this python and generate a plot for each cycle/member,
  and then look at the results, one by one, to double-check everything is ok.

DEPENDENCIES:
 See dependencies.py and the imports below.

AUTHOR and DATE:
 05/16/2022: Ricardo M. Campos, first version.

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

import matplotlib
matplotlib.use('Agg') # for backend plots, not for rendering in a window
import xarray as xr
import numpy as np
from pylab import *
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import sys
import pandas as pd
import cartopy.crs as ccrs
import cartopy
from matplotlib import ticker
import matplotlib.colors as colors
from matplotlib.colors import BoundaryNorm
# import pickle
import warnings
warnings.filterwarnings("ignore")
import wread
palette = plt.cm.jet
dpalette = plt.cm.RdBu_r

ftime=str(sys.argv[1]) # ftime="20000108" 
ftag=str(sys.argv[2]) # ftag="stream1_c00"

# lowest period (upper limit frequency) for the directional wave spectra (2D) polat plot
lper=4.5
# Paths
wrun="/work/noaa/marine/ricardo.campos/work/ww3runs/results/"+str(ftag.split('_')[0])+"/"+str(ftag.split('_')[1])
pbuoys="/work/noaa/marine/ricardo.campos/data/buoys/NDBC/ncformat/wparam"
opath="/work/noaa/marine/ricardo.campos/work/ww3runs/results/check/figs"
print(" ")

#  WW3 Fields ----------
fname=wrun+"/ww3gefs."+ftime+"_field.grib2"
ds = xr.open_dataset(fname, engine='cfgrib')
# ds.br.load_incrementally(time=10, nproc=4)
wtime1 = np.array(ds.time.values + ds.step.values )
lat = np.array(ds.latitude.values); lon = np.array(ds.longitude.values)
whs = np.array(ds['swh'].values[:,:,:]); units_whs = str(ds['swh'].units)
wuwnd = np.array(ds['u'].values[:,:,:]); units_wuwnd = str(ds['u'].units)
wvwnd = np.array(ds['v'].values[:,:,:]); units_wuwnd = str(ds['v'].units)
wtm = np.array(ds['mwp'].values[:,:,:]); units_wtm = str(ds['mwp'].units)
ds.close(); del ds, fname
# ----------

#   TIME-SERIES
# Observations NDBC ----------------
yrange=np.unique(np.array([np.int(np.str(wtime1[0])[0:4]),np.int(np.str(wtime1[-1])[0:4])]))
i=0; c=0
while c==0:

	stname=repr(np.int(41000)+np.int(i))
	try:
		#fname=str(pbuoys+"/"+repr(41000+i)+"h"+ftime[0:4]+".nc")
		#tso_time,lixo,tso_lat,tso_lon,lixo,lixo,lixo,lixo,lixo,tso_wnd,lixo,tso_hs,tso_tm,lixo,lixo  = wread.tseriesnc_ndbc(fname)
		tso_time=[]; tso_wnd=[]; tso_hs=[]; tso_tm=[]
		for j in range(0,yrange.shape[0]):
			fname=str(pbuoys+"/"+stname+"h"+repr(yrange[j])+".nc")
			atso_time,lixo,tso_lat,tso_lon,lixo,lixo,lixo,lixo,lixo,atso_wnd,lixo,atso_hs,atso_tm,lixo,lixo  = wread.tseriesnc_ndbc(fname)
			if j==0:
				tso_time=atso_time
				tso_wnd=atso_wnd
				tso_hs=atso_hs
				tso_tm=atso_tm
			else:
				tso_time=np.append(tso_time,atso_time)
				tso_wnd=np.append(tso_wnd,atso_wnd)
				tso_hs=np.append(tso_hs,atso_hs)
				tso_tm=np.append(tso_tm,atso_tm)

		# WAVEWATCH III ----------------
		fname=str(wrun+"/ww3gefs."+ftime+"_tab.nc")
		wtime2,lixo,lixo,lixo,tsm_hs,tsm_tm,lixo,tsm_dm,lixo,lixo,lixo = wread.tseriesnc_ww3(fname,stname)
	except:
		print(" "+stname+"h"+ftime[0:4]+".nc not available")
	else:
		ind = np.where( (tso_time>wtime1.min()) & (tso_time<wtime1.max()) & (tso_hs>0.1) )
		if size(ind)>50:
			c=1
			print(" ok "+stname)
			del fname

	i=i+1


#   SPECTRA
# WAVEWATCH III ----------------
fname=wrun+"/ww3gefs."+ftime+"_spec.nc"
wtime3,lixo,lixo,lixo,spm_freq,spm_freq1,spm_freq2,spm_dfreq,spm_pspec,spm_dmspec,spm_dire,spm_dspec,tsm_wnd,lixo = wread.spec_ww3(fname,stname)
if np.array_equal(wtime1,wtime3)==False: 
	sys.exit(' Error: _spec.nc with different time array from field.grib2.')
# ---------------------------------
print(" Read data OK")

# FIGURE ----------------

# for the 2D polar plot:
indf=int(np.where(abs(spm_freq-(1/lper))==min(abs(spm_freq-(1/lper))))[0][0])
ndire=np.zeros((spm_dire.shape[0]+2),'f'); ndire[1:-1]=spm_dire[:]; ndire[0]=0; ndire[-1]=360
angle = np.radians(ndire)
r, theta = np.meshgrid(spm_freq[0:indf], angle)
# ----------------------------------------

levelswnd = np.linspace(0,30,101)
levelshs = np.linspace(0,12,101)
levelstm = np.linspace(0,16,101)

fig, axs = plt.subplots(nrows=3,ncols=4,subplot_kw={'projection': ccrs.PlateCarree()},figsize=(19,10))
#   Fields -----------------
# First time step
gl = axs[0,0].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='grey', alpha=0.5, linestyle='--')
gl.xlabel_style = {'size': 9, 'color': 'k','rotation':0}; gl.ylabel_style = {'size': 9, 'color': 'k','rotation':0}
gl.ylabels_right = False; gl.xlabels_top = False
norm = BoundaryNorm(levelswnd, ncolors=palette.N, clip=False)
im = axs[0,0].pcolormesh(lon, lat, np.array(np.sqrt(wuwnd[0,:,:]**2 + wvwnd[0,:,:]**2)),shading='flat',cmap=palette,norm=norm, zorder=2)
axs[0,0].text(tso_lon, tso_lat,'*',color='k', size=12, zorder=3)
axs[0,0].add_feature(cartopy.feature.OCEAN,facecolor=("white"))
axs[0,0].add_feature(cartopy.feature.LAND,facecolor=("lightgrey"), edgecolor='grey',linewidth=0.5)
axs[0,0].add_feature(cartopy.feature.BORDERS, edgecolor='grey', linestyle='-',linewidth=0.5, alpha=1)
axs[0,0].coastlines(resolution='110m', color='grey',linewidth=0.5, linestyle='-', alpha=1)
axs[0,0].set_title(" 10mWind ("+units_wuwnd+") "+pd.to_datetime(wtime1[:][0]).strftime('%Y/%m/%d %H')+'Z') 
cax = axs[0,0].inset_axes([1.04, 0.2, 0.05, 0.6], transform=axs[0,0].transAxes)
cbar = plt.colorbar(im, ax=axs[0,0], cax=cax, extend='max')
tick_locator = ticker.MaxNLocator(nbins=6); cbar.locator = tick_locator; cbar.update_ticks()
#
gl = axs[0,1].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='grey', alpha=0.5, linestyle='--')
gl.xlabel_style = {'size': 9, 'color': 'k','rotation':0}; gl.ylabel_style = {'size': 9, 'color': 'k','rotation':0}
gl.ylabels_right = False; gl.xlabels_top = False
norm = BoundaryNorm(levelshs, ncolors=palette.N, clip=False)
im = axs[0,1].pcolormesh(lon, lat, whs[0,:,:],shading='flat',cmap=palette,norm=norm, zorder=2)
axs[0,1].text(tso_lon, tso_lat,'*',color='k', size=12, zorder=3)
axs[0,1].add_feature(cartopy.feature.OCEAN,facecolor=("white"))
axs[0,1].add_feature(cartopy.feature.LAND,facecolor=("lightgrey"), edgecolor='grey',linewidth=0.5)
axs[0,1].add_feature(cartopy.feature.BORDERS, edgecolor='grey', linestyle='-',linewidth=0.5, alpha=1)
axs[0,1].coastlines(resolution='110m', color='grey',linewidth=0.5, linestyle='-', alpha=1)
axs[0,1].set_title(" Hs ("+units_whs+") "+pd.to_datetime(wtime1[:][0]).strftime('%Y/%m/%d %H')+'Z') 
cax = axs[0,1].inset_axes([1.04, 0.2, 0.05, 0.6], transform=axs[0,1].transAxes)
cbar = plt.colorbar(im, ax=axs[0,1], cax=cax, extend='max')
tick_locator = ticker.MaxNLocator(nbins=6); cbar.locator = tick_locator; cbar.update_ticks()
#
gl = axs[0,2].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='grey', alpha=0.5, linestyle='--')
gl.xlabel_style = {'size': 9, 'color': 'k','rotation':0}; gl.ylabel_style = {'size': 9, 'color': 'k','rotation':0}
gl.ylabels_right = False; gl.xlabels_top = False
norm = BoundaryNorm(levelstm, ncolors=palette.N, clip=False)
im = axs[0,2].pcolormesh(lon, lat, wtm[0,:,:],shading='flat',cmap=palette,norm=norm, zorder=2)
axs[0,2].text(tso_lon, tso_lat,'*',color='k', size=12, zorder=3)
axs[0,2].add_feature(cartopy.feature.OCEAN,facecolor=("white"))
axs[0,2].add_feature(cartopy.feature.LAND,facecolor=("lightgrey"), edgecolor='grey',linewidth=0.5)
axs[0,2].add_feature(cartopy.feature.BORDERS, edgecolor='grey', linestyle='-',linewidth=0.5, alpha=1)
axs[0,2].coastlines(resolution='110m', color='grey',linewidth=0.5, linestyle='-', alpha=1)
axs[0,2].set_title(" Tm ("+units_wtm+") "+pd.to_datetime(wtime1[:][0]).strftime('%Y/%m/%d %H')+'Z')  
cax = axs[0,2].inset_axes([1.04, 0.2, 0.05, 0.6], transform=axs[0,2].transAxes)
cbar = plt.colorbar(im, ax=axs[0,2], cax=cax, extend='max')
tick_locator = ticker.MaxNLocator(nbins=6); cbar.locator = tick_locator; cbar.update_ticks()
# Last time step ----
gl = axs[1,0].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='grey', alpha=0.5, linestyle='--')
gl.xlabel_style = {'size': 9, 'color': 'k','rotation':0}; gl.ylabel_style = {'size': 9, 'color': 'k','rotation':0}
gl.ylabels_right = False; gl.xlabels_top = False
norm = BoundaryNorm(levelswnd, ncolors=palette.N, clip=False)
im = axs[1,0].pcolormesh(lon, lat, np.array(np.sqrt(wuwnd[-1,:,:]**2 + wvwnd[-1,:,:]**2)),shading='flat',cmap=palette,norm=norm, zorder=2)
axs[1,0].text(tso_lon, tso_lat,'*',color='k', size=12, zorder=3)
axs[1,0].add_feature(cartopy.feature.OCEAN,facecolor=("white"))
axs[1,0].add_feature(cartopy.feature.LAND,facecolor=("lightgrey"), edgecolor='grey',linewidth=0.5)
axs[1,0].add_feature(cartopy.feature.BORDERS, edgecolor='grey', linestyle='-',linewidth=0.5, alpha=1)
axs[1,0].coastlines(resolution='110m', color='grey',linewidth=0.5, linestyle='-', alpha=1)
axs[1,0].set_title(" 10mWind ("+units_wuwnd+") "+pd.to_datetime(wtime1[:][-1]).strftime('%Y/%m/%d %H')+'Z') 
cax = axs[1,0].inset_axes([1.04, 0.2, 0.05, 0.6], transform=axs[1,0].transAxes)
cbar = plt.colorbar(im, ax=axs[1,0], cax=cax, extend='max')
tick_locator = ticker.MaxNLocator(nbins=6); cbar.locator = tick_locator; cbar.update_ticks()
#
gl = axs[1,1].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='grey', alpha=0.5, linestyle='--')
gl.xlabel_style = {'size': 9, 'color': 'k','rotation':0}; gl.ylabel_style = {'size': 9, 'color': 'k','rotation':0}
gl.ylabels_right = False; gl.xlabels_top = False
norm = BoundaryNorm(levelshs, ncolors=palette.N, clip=False)
im = axs[1,1].pcolormesh(lon, lat, whs[-1,:,:],shading='flat',cmap=palette,norm=norm, zorder=2)
axs[1,1].text(tso_lon, tso_lat,'*',color='k', size=12, zorder=3)
axs[1,1].add_feature(cartopy.feature.OCEAN,facecolor=("white"))
axs[1,1].add_feature(cartopy.feature.LAND,facecolor=("lightgrey"), edgecolor='grey',linewidth=0.5)
axs[1,1].add_feature(cartopy.feature.BORDERS, edgecolor='grey', linestyle='-',linewidth=0.5, alpha=1)
axs[1,1].coastlines(resolution='110m', color='grey',linewidth=0.5, linestyle='-', alpha=1)
axs[1,1].set_title(" Hs ("+units_whs+") "+pd.to_datetime(wtime1[:][-1]).strftime('%Y/%m/%d %H')+'Z') 
cax = axs[1,1].inset_axes([1.04, 0.2, 0.05, 0.6], transform=axs[1,1].transAxes)
cbar = plt.colorbar(im, ax=axs[1,1], cax=cax, extend='max')
tick_locator = ticker.MaxNLocator(nbins=6); cbar.locator = tick_locator; cbar.update_ticks()
#
gl = axs[1,2].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='grey', alpha=0.5, linestyle='--')
gl.xlabel_style = {'size': 9, 'color': 'k','rotation':0}; gl.ylabel_style = {'size': 9, 'color': 'k','rotation':0}
gl.ylabels_right = False; gl.xlabels_top = False
norm = BoundaryNorm(levelstm, ncolors=palette.N, clip=False)
im = axs[1,2].pcolormesh(lon, lat, wtm[-1,:,:],shading='flat',cmap=palette,norm=norm, zorder=2)
axs[1,2].text(tso_lon, tso_lat,'*',color='k', size=12, zorder=3)
axs[1,2].add_feature(cartopy.feature.OCEAN,facecolor=("white"))
axs[1,2].add_feature(cartopy.feature.LAND,facecolor=("lightgrey"), edgecolor='grey',linewidth=0.5)
axs[1,2].add_feature(cartopy.feature.BORDERS, edgecolor='grey', linestyle='-',linewidth=0.5, alpha=1)
axs[1,2].coastlines(resolution='110m', color='grey',linewidth=0.5, linestyle='-', alpha=1)
axs[1,2].set_title(" Tm ("+units_wtm+") "+pd.to_datetime(wtime1[:][-1]).strftime('%Y/%m/%d %H')+'Z')  
cax = axs[1,2].inset_axes([1.04, 0.2, 0.05, 0.6], transform=axs[1,2].transAxes)
cbar = plt.colorbar(im, ax=axs[1,2], cax=cax, extend='max')
tick_locator = ticker.MaxNLocator(nbins=6); cbar.locator = tick_locator; cbar.update_ticks()
# ----
# Info
axs[2,3].remove()
fig.text(0.40,0.60,"Hs_mean : "+str(np.round(np.nanmean(whs[:,:,:]),2)),color='k', size=15, horizontalalignment='left', verticalalignment='center', transform=axs[2,3].transAxes)
fig.text(0.40,0.45,"Hs_max : "+str(np.round(np.nanmax(whs[:,:,:]),2)),color='k', size=15, horizontalalignment='left', verticalalignment='center', transform=axs[2,3].transAxes)
fig.text(0.40,0.30,"Tm_mean : "+str(np.round(np.nanmean(wtm[:,:,:]),2)),color='k', size=15, horizontalalignment='left', verticalalignment='center', transform=axs[2,3].transAxes)
fig.text(0.40,0.15,"Tm_max : "+str(np.round(np.nanmax(wtm[:,:,:]),2)),color='k', size=15, horizontalalignment='left', verticalalignment='center', transform=axs[2,3].transAxes)
fig.text(0.40,0.0,"Time Steps : "+str(np.int(wtime1.shape[0])),color='k', size=15, horizontalalignment='left', verticalalignment='center', transform=axs[2,3].transAxes)
# -----------------
# WW3 Wave Spectra -----------------
slevels = np.linspace(0.1,np.nanpercentile(spm_dspec,99.9),201)
#
ndspec=np.zeros((spm_freq.shape[0],ndire.shape[0]),'f')
ndspec[:,1:-1]=spm_dspec[0,:,:]
for i in range(0,spm_freq.shape[0]):
	ndspec[i,-1]=float((ndspec[i,-2]+ndspec[i,1])/2.)
	ndspec[i,0]=float((ndspec[i,-2]+ndspec[i,1])/2.)

axs[0,3].remove()
axs[0,3] = fig.add_subplot(3, 4, 4, projection='polar')
axs[0,3].set_theta_zero_location('N')
axs[0,3].set_theta_direction(-1)
axs[0,3].set_rlabel_position(-135)
axs[0,3].set_rticks([0.1,0.15,0.20]); axs[0,3].set_rmax(1/lper)
im = axs[0,3].contourf(theta, r, ndspec[0:indf,:].T,slevels,cmap=plt.cm.gist_stern_r,norm=colors.PowerNorm(gamma=0.5), extend="max")
axs[0,3].set_title('WW3Spec '+stname+', '+pd.to_datetime(wtime3[:][0]).strftime('%Y/%m/%d %H')+'Z',size=11)
cax = axs[0,3].inset_axes([1.13, 0.2, 0.05, 0.6], transform=axs[0,3].transAxes)
cbar = plt.colorbar(im, ax=axs[0,3], cax=cax)
tick_locator = ticker.MaxNLocator(nbins=6); cbar.locator = tick_locator; cbar.update_ticks()
del im, ndspec

ndspec=np.zeros((spm_freq.shape[0],ndire.shape[0]),'f')
ndspec[:,1:-1]=spm_dspec[-1,:,:]
for i in range(0,spm_freq.shape[0]):
	ndspec[i,-1]=float((ndspec[i,-2]+ndspec[i,1])/2.)
	ndspec[i,0]=float((ndspec[i,-2]+ndspec[i,1])/2.)

axs[1,3].remove()
axs[1,3] = fig.add_subplot(3, 4, 8, projection='polar')
axs[1,3].set_theta_zero_location('N')
axs[1,3].set_theta_direction(-1)
axs[1,3].set_rlabel_position(-135)
axs[1,3].set_rticks([0.1,0.15,0.20]); axs[1,3].set_rmax(1/lper)
im = axs[1,3].contourf(theta, r, ndspec[0:indf,:].T,slevels,cmap=plt.cm.gist_stern_r,norm=colors.PowerNorm(gamma=0.5), extend="max")
axs[1,3].set_title('WW3Spec '+stname+', '+pd.to_datetime(wtime3[:][-1]).strftime('%Y/%m/%d %H')+'Z',size=11)
cax = axs[1,3].inset_axes([1.13, 0.2, 0.05, 0.6], transform=axs[1,3].transAxes)
cbar = plt.colorbar(im, ax=axs[1,3], cax=cax)
tick_locator = ticker.MaxNLocator(nbins=6); cbar.locator = tick_locator; cbar.update_ticks()
del im, ndspec
# -----------------
# Time-Series -----------------
axs[2,0].remove(); axs[2,1].remove(); axs[2,2].remove()
iaux=np.intersect1d(wtime3, np.array(tso_time, dtype='datetime64[h]'), assume_unique=False, return_indices=True)
axs[2,0] = fig.add_subplot(3, 4, 9, projection='rectilinear')
axs[2,0].plot_date(wtime3,tsm_wnd,color='b', linestyle='-',marker='',linewidth=2.0, label="GEFS", zorder=3)
if size(iaux)>0:
	axs[2,0].plot_date(tso_time[np.min(iaux[2]):np.max(iaux[2])+1],tso_wnd[np.min(iaux[2]):np.max(iaux[2])+1],'k.',label='buoy', zorder=2)

axs[2,0].xaxis.set_major_formatter( DateFormatter('%d') ); axs[2,0].fmt_xdata = DateFormatter('%d')
axs[2,0].grid(linewidth=0.5, color='grey', alpha=0.5, linestyle='--', zorder=1)
axs[2,0].legend(loc='best', fontsize=9)
axs[2,0].set_xlabel('Date', fontsize=12); axs[2,0].set_ylabel(" 10mWind ("+units_wuwnd+") ", fontsize=12) 
axs[2,0].axis('tight'); del iaux
#
iaux=np.intersect1d(wtime2, np.array(tso_time, dtype='datetime64[h]'), assume_unique=False, return_indices=True)
axs[2,1] = fig.add_subplot(3, 4, 10, projection='rectilinear')
axs[2,1].plot_date(wtime2,tsm_hs,color='b', linestyle='-',marker='',linewidth=2.0, label="WW3GEFS", zorder=3)
if size(iaux)>0:
	axs[2,1].plot_date(tso_time[np.min(iaux[2]):np.max(iaux[2])+1],tso_hs[np.min(iaux[2]):np.max(iaux[2])+1],'k.',label='buoy', zorder=2)

axs[2,1].xaxis.set_major_formatter( DateFormatter('%d') ); axs[2,1].fmt_xdata = DateFormatter('%d')
axs[2,1].grid(linewidth=0.5, color='grey', alpha=0.5, linestyle='--', zorder=1)
axs[2,1].legend(loc='best', fontsize=9)
axs[2,1].set_xlabel('Date', fontsize=12); axs[2,1].set_ylabel(" Hs ("+units_whs+") ", fontsize=12) 
axs[2,1].axis('tight'); del iaux
#
iaux=np.intersect1d(wtime2, np.array(tso_time, dtype='datetime64[h]'), assume_unique=False, return_indices=True)
axs[2,2] = fig.add_subplot(3, 4, 11, projection='rectilinear')
axs[2,2].plot_date(wtime2,tsm_tm,color='b', linestyle='-',marker='',linewidth=2.0, label="WW3GEFS", zorder=3)
if size(iaux)>0:
	axs[2,2].plot_date(tso_time[np.min(iaux[2]):np.max(iaux[2])+1],tso_tm[np.min(iaux[2]):np.max(iaux[2])+1],'k.',label='buoy', zorder=2)

axs[2,2].xaxis.set_major_formatter( DateFormatter('%d') ); axs[2,2].fmt_xdata = DateFormatter('%d')
axs[2,2].grid(linewidth=0.5, color='grey', alpha=0.5, linestyle='--', zorder=1)
axs[2,2].legend(loc='best', fontsize=9)
axs[2,2].set_xlabel('Date', fontsize=12); axs[2,2].set_ylabel(" Tm ("+units_wtm+") ", fontsize=12) 
axs[2,2].axis('tight'); del iaux

fig.canvas.draw() # https://github.com/SciTools/cartopy/issues/1207
fig.tight_layout()
plt.savefig(opath+'/CheckWW3GEFS_'+ftime+'_'+str(ftag.split('_')[1])+'.png', dpi=200, facecolor='w', edgecolor='w',
		orientation='portrait', papertype=None, format='png',transparent=False, pad_inches=0.1)

plt.close('all'); del axs, fig
print(" Plot OK")

