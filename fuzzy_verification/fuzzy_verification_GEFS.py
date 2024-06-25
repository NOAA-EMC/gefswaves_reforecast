#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
fuzzy_verification_GEFS.py

VERSION AND LAST UPDATE:
 v1.0  10/30/2023
 v1.1  11/03/2023
 v1.2  05/20/2024

PURPOSE:
 Validation of long-term probabilistic wave forecasts using 
 fuzzy verification. This can be used to assess and optimize
 the parameters of operational probability maps.
 The probability maps are based on GEFSv12 and have been validated against
  NDBC buoys and GEFS hindcasts by appending 24-hour segments 
  from consecutive cycles.

USAGE:
 Four input arguments are required:
 - Station ID (41048 etc)
 - Forecast Lead Time (Day), Initial (ex. 7)
 - Forecast Lead Time (Day), Final (ex. 14)
 - Output path where output files will be saved.

DEPENDENCIES:
 See the imports below.

AUTHOR and DATE:
 10/30/2023: Ricardo M. Campos, first version, named fuzzy_verification_ProbMaps.py
 11/03/2023: Ricardo M. Campos, fuzzy verification including categorical 
  and probabilistic assessments.
 05/20/2024: Ricardo M. Campos, including all the station in the analysis,
  GDAS has been excluded, and code renamed to fuzzy_verification_GEFS.py

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

import matplotlib
matplotlib.use('Agg')
import numpy as np
import sys
import pandas as pd
import matplotlib.pyplot as plt
import yaml
from pvalstats import ModelObsPlot
from wcpval import *
import warnings; warnings.filterwarnings("ignore")

sl=13
matplotlib.rcParams.update({'font.size': sl}); plt.rc('font', size=sl) 
matplotlib.rc('xtick', labelsize=sl); matplotlib.rc('ytick', labelsize=sl); matplotlib.rcParams.update({'font.size': sl})

if __name__ == "__main__":

    # point manual selection (index)
    # All
    bid=np.arange(0,21).astype('int')
    # Gulf
    bid=np.array([8,9,10,11]).astype('int')
    # Atlantic
    bid=np.array([0,1,2,3,4,5,6,7,12]).astype('int')
    # Pacific
    bid=np.array([13,14,15,16,17,18,19,20]).astype('int')
    # Tropical
    bid=np.array([2,3,4,5,6,7,8,9,10,11]).astype('int')
    # Extratropical
    bid=np.array([12,13,14,15,16,17,18,19]).astype('int')

    lbid = len(bid)
    # variable (u10 or hs)
    wvar='hs'
    # Forecast Lead Time (Day) and intervall
    ltime1=7
    ltime2=14
    # output path
    opath="/home/ricardo/cimas/analysis/3assessments/fuzzy_verification/output"
    # opath="/work/noaa/marine/ricardo.campos/work/analysis/3assessments/fuzzy_verification/output"
    # file tag to same output fils
    ftag=opath+"/Validation_"+wvar+"_"
    # ---- Read statistical parameters (exaclty the same as the operational config file) -----
    print(" "); print(" Reading yaml configuration file ...")
    with open('probmaps_gefs.yaml', 'r') as file:
        wconfig = yaml.safe_load(file)

    if wvar=="u10":
        fqlev = np.array(wconfig['qlev_wnd']).astype('float')
        qlev = fqlev/1.94 # m/s
    else:
        qlev = np.array(wconfig["qlev_"+wvar]).astype('float')
        fqlev = qlev

    plevels = np.array(wconfig['hplevels']).astype('float'); lplev=int(len(plevels)-1)
    nmax = int(wconfig['nmax']); spws = int(wconfig['spws']); spctl = float(wconfig['spctl'])

    # -------------------------------------

    # READ DATA
    # list of netcdf files generated with buildfuzzydataset.py (GEFS, GDAS, and NDBC buoy)
    # ls -d $PWD/*.nc > list.txt &
    wlist = np.atleast_1d(np.loadtxt('list.txt',dtype=str)) 
    gdata = read_data(wlist,bid,ltime1,ltime2,wvar)
    indlat = gdata['indlat']; indlon = gdata['indlon']

    # lengths
    lstw=int(len(wlist))
    lft = gdata['lft']
    lct = len(gdata['ctime'])
    lensm = len(gdata['ensm'])
    llatm = len(gdata['latm'])
    llonm = len(gdata['lonm'])

    # Reshape
    # join bstations and forecast time
    ndbc = np.array(gdata['ndbc'].reshape(lct,lbid*lft))
    gefs_hindcast = gdata['gefs_hindcast'].reshape(lct,lbid*lft,lensm,llatm,llonm)
    gefs_forecast = gdata['gefs_forecast'].reshape(lct,lbid*lft,lensm,llatm,llonm)
    # join cycles
    r_ndbc = np.array(ndbc.reshape(lct*lbid*lft))
    r_gefs_hindcast = gefs_hindcast.reshape((lct*lbid*lft,lensm,llatm,llonm))
    r_gefs_forecast = gefs_forecast.reshape((lct*lbid*lft,lensm,llatm,llonm))

    # ===== Mimic the operational probability maps methodology ===== 

    # Ground truth (nmax within each cycle), buoy and Hindcast:
    ndbc_t = nmaxsel(ndbc,nmax)
    # model
    gefs_hindcast_t = np.zeros((lstw),'f')*np.nan
    for i in range(0,lstw):
        aux = np.nanmean(gefs_hindcast[i,:,:,indlat,indlon],axis=1) # average ensemble members (hindcast)
        ind = np.where(aux>-999.)
        if np.size(ind)>int(np.floor(lft/2)) :
            ind = ind[0]
            gefs_hindcast_t[i] = np.nanmean(np.sort(aux[ind])[-nmax::])
        else:
            print(" forecast time series incomplete "+repr(i))
            gefs_hindcast_t[i] = np.nan

        del ind,aux


    #  -- Probabilistic Forecast array --  
    gspws=int(np.floor(spws/np.diff(gdata['latm']).mean())/2)
    prob_gefs_forecast, fmod_result = probforecast(nmax,gspws,spctl,gefs_forecast,qlev)


    #  ++++++++ VALIDATION +++++++++++++++

    # Quality of ground truth: Initial Validation Hindcast X Buoy, deterministic
    for i in range(0,len(qlev)):
        # Deterministic, matching positions and time
        if i==0:
            mop=ModelObsPlot(np.nanmean(r_gefs_hindcast[:,:,indlat,indlon],axis=1),r_ndbc,ftag=ftag+"DeterministicVal_GEFShindcast_")
            mop.scatterplot(); mop.pdf(); mop.qqplot()

    # ===== Categorical approach, verification of binary events =====

    # Binary Categorical array
    fprob_gefs_forecast = probforecast_binary(prob_gefs_forecast,qlev,plevels)
    # Validation scores
    for i in range(0,len(qlev)):
        for j in range(0,lplev):
            ceval_gefs_ndbc = categorical_bin_eval(fprob_gefs_forecast[:,i,j],ndbc_t,qlev[i])
            ceval_gefs_ghnd = categorical_bin_eval(fprob_gefs_forecast[:,i,j],gefs_hindcast_t,qlev[i])
            dict_list = [
                {'Name': 'gefs_ndbc', **ceval_gefs_ndbc},
                {'Name': 'gefs_gefsHindcast', **ceval_gefs_ghnd}
            ]

            # Convert the list of dictionaries into a pandas DataFrame
            df = pd.DataFrame(dict_list)
            df.to_csv(ftag+"CategoricalVal_Lev"+repr(np.round(fqlev[i],2))+"_Prob"+repr(np.round(plevels[j],2))+".csv",sep='\t', index=False)


    # Probabilistic approach, verification of event frequency
    # BS, CRPS, ROC, reliability curve

    onames = ['gefs_ndbc', 'gefs_gefsHindcast']

    plt.close('all')
    # --- Brier Score ---
    for i in range(0,len(qlev)):
        fbriers=[]
        fftag=ftag+"ProbEvents_NDBC_"+repr(np.round(fqlev[i],2))
        briers = brier_score(prob_gefs_forecast[:,i],ndbc_t,qlev[i],gdata['cdate'],fftag)
        fbriers = np.append(fbriers,briers)

        fftag=ftag+"ProbEvents_GEFShindcast_"+repr(np.round(fqlev[i],2))
        briers = brier_score(prob_gefs_forecast[:,i],gefs_hindcast_t,qlev[i],gdata['cdate'],fftag)
        fbriers = np.append(fbriers,briers)

        bdata = {'Name': onames, 'Result': np.round(fbriers,3)}
        df = pd.DataFrame(bdata)
        df.to_csv(ftag+"ProbabilisticVal_BrierScore_Lev"+repr(np.round(fqlev[i],2))+".csv",sep='\t', index=False)

    plt.close('all')

    # --- CRPS ---
    fcrps_mean=[]
    fftag=ftag+"ProbEvents_NDBC"
    mcrps, mcrps_mean = crps(fmod_result,ndbc_t,gdata['cdate'],fftag)
    fcrps_mean = np.append(fcrps_mean,mcrps_mean)

    fftag=ftag+"ProbEvents_GEFShindcast"
    mcrps, mcrps_mean  = crps(fmod_result,gefs_hindcast_t,gdata['cdate'],fftag)
    fcrps_mean = np.append(fcrps_mean,mcrps_mean)

    bdata = {'Name': onames, 'Result': np.round(fcrps_mean,4)}
    df = pd.DataFrame(bdata)
    df.to_csv(ftag+"ProbabilisticVal_CRPS.csv",sep='\t', index=False)

    plt.close('all')

    # --- ROC Curve ---
    for i in range(0,len(qlev)):
        froc=[]
        fftag=ftag+"ProbEvents_NDBC_"+repr(np.round(fqlev[i],2))
        true_binary = (ndbc_t > qlev[i]).astype(int)
        roc = roc_plot(true_binary,prob_gefs_forecast[:,i],fftag)
        froc = np.append(froc,roc); del true_binary, fftag

        fftag=ftag+"ProbEvents_GEFShindcast_"+repr(np.round(fqlev[i],2))
        true_binary = (gefs_hindcast_t > qlev[i]).astype(int)
        roc = roc_plot(true_binary,prob_gefs_forecast[:,i],fftag)
        froc = np.append(froc,roc); del true_binary, fftag

        bdata = {'Name': onames, 'Result': np.round(froc,3)}
        df = pd.DataFrame(bdata)
        df.to_csv(ftag+"ProbabilisticVal_ROCauc_Lev"+repr(np.round(fqlev[i],2))+".csv",sep='\t', index=False)

    plt.close('all')

    # --- Reliability Curve ---
    for i in range(0,len(qlev)):

        fftag=ftag+"ProbEvents_NDBC_"+repr(np.round(fqlev[i],2))
        true_binary = (ndbc_t > qlev[i]).astype(int)
        reliability_curve(true_binary,prob_gefs_forecast[:,i],fftag)

        fftag=ftag+"ProbEvents_GEFShindcast_"+repr(np.round(fqlev[i],2))
        true_binary = (gefs_hindcast_t > qlev[i]).astype(int)
        reliability_curve(true_binary,prob_gefs_forecast[:,i],fftag)

    plt.close('all')

