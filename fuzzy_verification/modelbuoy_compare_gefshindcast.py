#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
modelbuoy_compare_gefshindcast.py

VERSION AND LAST UPDATE:
 v1.0  06/24/2025

PURPOSE:
 Comparison between GEFS (hindcasted), after extracting points with extractPoint_makeHindcast_GEFS.py, and buoys.

USAGE:
 Requires NDBC buoy data previously downloaded, and netcdf file with GEFS (output of extractPoint_makeHindcast_GEFS.py)
 Comparisons are for one month, selected.

OUTPUT:
 Png plots of Hs and Tp.

DEPENDENCIES:
 See setup.py and the imports below.

AUTHOR and DATE:
 06/24/2025: Ricardo M. Campos, first version.

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

import matplotlib
matplotlib.use('Agg')
import warnings; warnings.filterwarnings("ignore")
import numpy as np
from matplotlib.mlab import *
from pylab import *
import xarray as xr
import netCDF4 as nc
import time
from time import strptime
from calendar import timegm
import wread
import pvalstats
import mvalstats
from pvalstats import ModelObsPlot
# netcdf format
fnetcdf="NETCDF4"

sl=13
matplotlib.rcParams.update({'font.size': sl}); plt.rc('font', size=sl)
matplotlib.rc('xtick', labelsize=sl); matplotlib.rc('ytick', labelsize=sl); matplotlib.rcParams.update({'font.size': sl})


def timeseries(model,buoy,mdate,wlabel,ftag):
    '''
    Time-series plots of model (lines) and observations (discrete points).
    The x-axis is time
    '''

    merr=mvalstats.metrics(model[0,:],buoy)
    year=str(mdate[0].astype('datetime64[Y]').astype(int) + 1970)

    fig1, ax = plt.subplots(figsize=(9, 4))
    ax.plot(mdate, buoy, color='k', marker='.', linestyle='', linewidth=2., label='Obs', zorder=3)
    ax.plot(mdate, model[0,:], color='blue', linestyle='-', linewidth=2., label='Model', zorder=2)

    ax.fill_between(mdate,np.nanmin(model,axis=0), np.nanmax(model,axis=0), color='lightblue', alpha=0.6, zorder=1)

    ax.set_xlim(mdate[0], mdate[-1])
    ax.xaxis.set_major_formatter(DateFormatter('%b%d'))
    ax.fmt_xdata = DateFormatter('%b%d')

    ax.text(0.01, 0.99, "Bias "+str(np.round(merr[0],3)), transform=ax.transAxes, fontsize=12, verticalalignment='top', horizontalalignment='left')

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(sl - 2)

    ax.set_xlabel("Time ("+year+")")
    ax.set_ylabel(wlabel)
    ax.legend(fontsize=sl - 3)

    # plt.tight_layout()
    plt.grid(c='grey', ls='--', alpha=0.3, zorder=1)
    plt.savefig("TimeSeries_"+ftag+".png", dpi=200, facecolor='w', edgecolor='w', orientation='portrait',
        format='png', bbox_inches='tight', pad_inches=0.1)

    plt.close(fig1)


if __name__ == "__main__":

    # Paths
    ndbcp="/work/noaa/marine/ricardo.campos/data/buoys/NDBC/ncformat/wparam"
    # ndbcp="/home/ricardo/work/noaa/analysis/Week2ProbForecast/3assessments/fuzzy_verification/events/Pacific"

    fmodname="/work/noaa/marine/ricardo.campos/work/analysis/3assessments/events/GEFS.PointExtract.20240101to20240201.nc"
    # fmodname="/home/ricardo/work/noaa/analysis/Week2ProbForecast/3assessments/fuzzy_verification/events/Pacific/GEFS.PointExtract.20240101to20240201.nc"

    f=nc.Dataset(fmodname); ds = xr.open_dataset(fmodname)
    mtime=f.variables['time'][:]; mdate = ds.time.values[:]
    stname=f.variables['buoyID'][:]
    mhs=f.variables['gefs_ww3'][0,:,:,:] # HTSGW_surface
    mtp=f.variables['gefs_ww3'][1,:,:,:] # PERPW_surface
    f.close(); del f, ds

    bhs = np.zeros((len(stname),len(mtime)),'f')*np.nan
    btp = np.zeros((len(stname),len(mtime)),'f')*np.nan
    year=str(mdate[0].astype('datetime64[Y]').astype(int) + 1970)

    for i in range(0,len(stname)):

        try:
            res = wread.tseriesnc_ndbc(fname=ndbcp+"/"+stname[i]+"h"+year+".nc",anh=None)
            ahs = res['hs']
            atp = res['tp']
            btime = res['time']
            bdate = res['date']
            del res
        except:
            print(" Cannot open "+stname[i])
        else:

            if np.any(mhs[i,0,:]>0.):

                # First layer of simple quality-control (Obs)
                indq=np.where((ahs>30.)|(ahs<0.0))
                if np.size(indq)>0:
                    ahs[indq]=np.nan; del indq

                indq=np.where((atp>40.)|(atp<0.0))
                if np.size(indq)>0:
                    atp[indq]=np.nan; del indq

                c=0
                for t in range(0,len(mtime)):
                    indt=np.where(np.abs(btime-mtime[t])<1800.)
                    if np.size(indt)>0:
                        bhs[i,t] = float(np.nanmean(ahs[indt[0]]))
                        btp[i,t] = float(np.nanmean(atp[indt[0]]))
                        c=c+1
                        del indt

                if np.any(bhs[i,:]>0.):
                    # Hs
                    ftag="Hs_"+stname[i]+"_"
                    mop=ModelObsPlot(mhs[i,0,:],bhs[i,:],ftag=ftag)
                    mop.qqplot()
                    wlabel="Hs (m)"
                    timeseries(mhs[i,:,:],bhs[i,:],mdate,wlabel,ftag)
                
                    # Tp
                    ftag="Tp_"+stname[i]+"_"
                    mop=ModelObsPlot(mtp[i,0,:],btp[i,:],ftag=ftag)
                    mop.qqplot()
                    wlabel="Tp (s)"
                    timeseries(mtp[i,:,:],btp[i,:],mdate,wlabel,ftag)

                    print(" Ok "+stname[i])

