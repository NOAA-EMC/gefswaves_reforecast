#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plot_bpositions.py

VERSION AND LAST UPDATE:
 v1.0  05/15/2024

PURPOSE:
 Small script to plot the position of the buoys and the model grid points
  selected.

"""

import matplotlib
# matplotlib.use('Agg')
import netCDF4 as nc
from matplotlib import ticker
from pylab import *
import cartopy
import cartopy.crs as ccrs
import warnings; warnings.filterwarnings("ignore")

sl=14
matplotlib.rcParams.update({'font.size': sl}); plt.rc('font', size=sl) 
matplotlib.rc('xtick', labelsize=sl); matplotlib.rc('ytick', labelsize=sl); matplotlib.rcParams.update({'font.size': sl})

if __name__ == "__main__":

    fname = '/media/ricardo/ssdrmc/analysis/3assessments/fuzzy_verification/data/GEFS.BUOY.PointExtract.20220101.nc'

    # ==== READ ====
    f=nc.Dataset(fname)
    bid=f.variables['buoyID'][:]
    blat=f.variables['buoy_lat'][:]; blon=f.variables['buoy_lon'][:]
    lat=f.variables['lat'][:]; lon=f.variables['lon'][:]
    f.close(); del f

    # ==== Plot point positions ====
    # Buoys
    plt.figure(figsize=(7,5))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=-90))
    ax.set_extent([-165,-45,0.,65.], crs=ccrs.PlateCarree())
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='grey', alpha=0.5, linestyle='--')
    gl.xlabel_style = {'size': 9, 'color': 'k','rotation':0}; gl.ylabel_style = {'size': 9, 'color': 'k','rotation':0}
    ax.add_feature(cartopy.feature.OCEAN,facecolor=("white"))
    ax.add_feature(cartopy.feature.LAND,facecolor=("lightgrey"), edgecolor='grey',linewidth=0.5, zorder=2)
    ax.add_feature(cartopy.feature.BORDERS, edgecolor='grey', linestyle='-',linewidth=0.5, alpha=1, zorder=3)
    ax.coastlines(resolution='50m', color='dimgrey',linewidth=0.5, linestyle='-', alpha=1, zorder=4)
    # Buoy
    for i in range(0,len(bid)):
        ax.scatter(blon[i],blat[i],s=5, c='k', marker='o',transform=ccrs.PlateCarree(),zorder=3)
        for j in range(0,lon.shape[1]):
            for k in range(0,lon.shape[1]):
                ax.scatter(lon[i,j],lat[i,k],s=1, c='lightblue',alpha=0.5, marker='.',transform=ccrs.PlateCarree(),zorder=1)

    plt.savefig('Buoy_Positions.png', dpi=300, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
    plt.close()

