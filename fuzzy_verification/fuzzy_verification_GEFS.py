#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
fuzzy_verification_GEFS.py

VERSION AND LAST UPDATE:
 v1.0  10/30/2023
 v1.1  11/03/2023
 v1.2  05/20/2024
 v1.3  01/10/2025

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
 01/10/2025: Ricardo M. Campos, climatology and persistence added to the prob validation. Skill 
  scores BSS and CRPSS included. Unnecessary plots have been removed."

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
    # bid=np.arange(0,22).astype('int')
    bid = np.array([7,8,9,10])
    lbid = len(bid)
    # variable (u10 or hs)
    wvar='hs'
    # Forecast Lead Time (Day) and intervall
    ltime1=7
    ltime2=14
    # output path
    opath="/home/ricardo/work/noaa/analysis/Week2ProbForecast/3assessments/fuzzy_verification/output/Pacific/bc/np_tropical"
    # opath="/work/noaa/marine/ricardo.campos/work/analysis/3assessments/fuzzy_verification/output"
    # file tag for output file names
    ftag=opath+"/Validation_"+wvar+"_"

    # Cluster data, altogether: yes (1) or no (0)
    cpoints = 0

    # Bias correction of ground truth
    qmm_slope = 1.05
    qmm_intercept = 0. 

    # ---- Read config file and statistical parameters (exaclty the same as the operational) -----
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
    print(" Read yaml configuration file, OK")

    # -------------------------------------
    # READ DATA
    print(" Reading Model Data ...")
    # list of netcdf files generated with buildfuzzydataset.py
    # ls -d $PWD/*.nc > list.txt &
    wlist = np.atleast_1d(np.loadtxt('list.txt',dtype=str)) 
    gdata = read_data(wlist,bid,ltime1,ltime2,wvar)
    indlat = gdata['indlat']; indlon = gdata['indlon']
    print(" Reading Model Data, OK")

    print(" Prepare arrays ...")
    # lengths
    lstw=int(len(wlist))
    lft = gdata['lft']
    lct = len(gdata['ctime'])
    lensm = len(gdata['ensm'])
    llatm = len(gdata['latm'])
    llonm = len(gdata['lonm'])

    # Reshape
    if cpoints == 0:
        gtloop=int(lct*lbid)
        gefs_hindcast = gdata['gefs_hindcast'].reshape(gtloop,lft,lensm)
        gefs_forecast = gdata['gefs_forecast'].reshape(gtloop,lft,lensm,llatm,llonm)
        gefs_forecast_p = gdata['gefs_forecast_p'].reshape(gtloop,lft,lensm,llatm,llonm)
    else:
        gtloop=lct
        gefs_hindcast = gdata['gefs_hindcast'].reshape(gtloop,lbid*lft,lensm)
        gefs_forecast = gdata['gefs_forecast'].reshape(gtloop,lbid*lft,lensm,llatm,llonm)
        gefs_forecast_p = gdata['gefs_forecast_p'].reshape(gtloop,lbid*lft,lensm,llatm,llonm)

    print(" Prepare arrays, OK")

    # ===== Mimic the operational probability maps methodology ===== 
    print(" Compute probabilities (same as operational) ...")

    # Ground truth 
    gefs_hindcast_t = np.zeros((gtloop),'f')*np.nan
    for i in range(0,gtloop):
        aux = np.nanmean(gefs_hindcast[i,:,:],axis=1) # average ensemble members (hindcast)
        ind = np.where(aux>-999.)
        if np.size(ind)>int(np.floor(lft/2)) :
            ind = ind[0]
            gefs_hindcast_t[i] = np.nanmean(np.sort(aux[ind])[-nmax::])
        else:
            print(" forecast time series incomplete "+repr(i))
            gefs_hindcast_t[i] = np.nan

        del ind,aux

    # Apply calibration function (from QMM bias correction)
    gefs_hindcast_t = gefs_hindcast_t*qmm_slope + qmm_intercept



    #  -- Probabilistic Forecast array --  
    gspws=int(np.floor(spws/np.diff(gdata['latm']).mean())/2)
    prob_gefs_forecast, fmod_result = probforecast(nmax,gspws,spctl,gefs_forecast,qlev)
    prob_gefs_forecast_p, fmod_result_p = probforecast(nmax,gspws,spctl,gefs_forecast_p,qlev)

    # Climatology
    prob_gefs_forecast_c = np.array(prob_gefs_forecast*0.)
    aux = np.array( gefs_forecast[:,:,:,indlat,indlon] ).reshape(gefs_forecast.shape[0]*gefs_forecast.shape[1]*gefs_forecast.shape[2])
    for i in range(0,len(qlev)):
        ind = np.where(aux>qlev[i])
        if np.size(ind)>0:
            prob_gefs_forecast_c[:,i] = float(len(ind[0])/len(aux))
            del ind

    del aux
    fmod_result_c = np.array(fmod_result*0.)+np.nanmean(fmod_result)
    print(" Compute probabilities (same as operational), OK")


    #  ++++++++ VALIDATION +++++++++++++++
    print(" VALIDATION ... ")

    # ===== Categorical approach, verification of binary events =====

    # ---- Binary Categorical array ------
    print(" Binary Categorical ... ")

    # GEFS
    fprob_gefs_forecast = probforecast_binary(prob_gefs_forecast,qlev,plevels)
    # Validation scores
    for i in range(0,len(qlev)):
        for j in range(0,lplev):
            ceval_gefs_ghnd = categorical_bin_eval(fprob_gefs_forecast[:,i,j],gefs_hindcast_t,qlev[i])
            dict_list = [ {'Name': 'gefs_gefsHindcast', **ceval_gefs_ghnd} ]

            # Convert the list of dictionaries into a pandas DataFrame
            df = pd.DataFrame(dict_list)
            df.to_csv(ftag+"CategoricalVal_Lev"+repr(np.round(fqlev[i],2))+"_Prob"+repr(np.round(plevels[j],2))+".csv",sep='\t', index=False)

    # GEFS persistence
    fprob_gefs_forecast_p = probforecast_binary(prob_gefs_forecast_p,qlev,plevels)
    # Validation scores
    for i in range(0,len(qlev)):
        for j in range(0,lplev):
            ceval_gefs_ghnd_p = categorical_bin_eval(fprob_gefs_forecast_p[:,i,j],gefs_hindcast_t,qlev[i])
            dict_list = [ {'Name': 'gefs_gefsHindcast', **ceval_gefs_ghnd_p} ]

            # Convert the list of dictionaries into a pandas DataFrame
            df = pd.DataFrame(dict_list)
            df.to_csv(ftag+"CategoricalVal_Persistence_Lev"+repr(np.round(fqlev[i],2))+"_Prob"+repr(np.round(plevels[j],2))+".csv",sep='\t', index=False)

    print(" Binary Categorical, OK")

    # ========================================================
    # Probabilistic approach, verification of event frequency
    # BSS, CRPSS, ROC, Reliability Curve

    print(" Probabilistic validation.")
    print(" Brier skill score ...")

    plt.close('all')

    # --- Brier Skill Score ---
    # GEFS
    for i in range(0,len(qlev)):
        try:
            fbriers=[]
            briers = brier_score(prob_gefs_forecast[:,i],gefs_hindcast_t,qlev[i])
            briers_p = brier_score(prob_gefs_forecast_p[:,i],gefs_hindcast_t,qlev[i])
            briers_c = brier_score(prob_gefs_forecast_c[:,i],gefs_hindcast_t,qlev[i])

            bss_p = - (briers - briers_p)/briers_p
            fbriers = np.append(fbriers,bss_p)
            bss_c = - (briers - briers_c)/briers_c
            fbriers = np.append(fbriers,bss_c)

            bdata = {'Name': ['BSS_P','BSS_C'], 'Result': np.round(fbriers,3)}
            df = pd.DataFrame(bdata)
            df.to_csv(ftag+"ProbabilisticVal_BSS_Lev"+repr(np.round(fqlev[i],2))+".csv",sep='\t', index=False)
        except:
		    print("  - Did not compute Brier skill score "+repr(qlev[i]))
        else:
            print("  - Brier skill score for "+repr(qlev[i])+" ok")

    print(" Brier skill score, OK")

    # --- CRPS ---
    print(" CRPS ...")

    fcrps=[]
    mcrps, mcrps_mean  = crps(fmod_result,gefs_hindcast_t)
    mcrps_p, mcrps_mean_p  = crps(fmod_result_p,gefs_hindcast_t)
    mcrps_c, mcrps_mean_c  = crps(fmod_result_c,gefs_hindcast_t)

    crpss_p = 1 - (mcrps_mean/mcrps_mean_p)
    fcrps = np.append(fcrps,crpss_p)
    crpss_c = 1 - (mcrps_mean/mcrps_mean_c)
    fcrps = np.append(fcrps,crpss_c)

    bdata = {'Name': ['CRPSS_P','CRPSS_C'], 'Result': np.round(fcrps,4)}
    df = pd.DataFrame(bdata)
    df.to_csv(ftag+"ProbabilisticVal_CRPSS.csv",sep='\t', index=False)

    print(" CRPS, OK")

    mlabels = ['GEFS', 'GEFS_p', 'GEFS_c']
    # --- ROC Curve ---
    print(" ROC curve ...")
    for i in range(0,len(qlev)-1):
        try:
            prob_results = np.array([prob_gefs_forecast[:,i],prob_gefs_forecast_p[:,i],prob_gefs_forecast_c[:,i]])

            fftag=ftag+"ProbEvents_GEFShindcast_"+repr(np.round(fqlev[i],2))
            true_binary = (gefs_hindcast_t > qlev[i]).astype(int)
            roc = roc_plot2(true_binary,prob_results,mlabels,plevels,fftag)

            bdata = {'Name': mlabels, 'Result': np.round(roc,3)}
            df = pd.DataFrame(bdata)
            df.to_csv(ftag+"ProbabilisticVal_GEFShindcast_ROCauc_Lev"+repr(np.round(fqlev[i],2))+".csv",sep='\t', index=False)
        except:
		    print("  - Did not generate ROC curve for "+repr(qlev[i]))
        else:
            print("  - ROC curve for "+repr(qlev[i])+" ok")

    plt.close('all')
    print(" ROC curve, OK")

    # --- Reliability Curve ---
    print(" Reliability Curve ...")
    for i in range(0,len(qlev)):
        try:

            if i<2:
                nbins=10
                pmax=1
            else:
                nbins=5
                pmax=0.75

            prob_results = np.array([prob_gefs_forecast[:,i],prob_gefs_forecast_p[:,i],prob_gefs_forecast_c[:,i]])

            ght = gefs_hindcast_t[~np.isnan(gefs_hindcast_t)]  
            ind=np.where(ght >=fqlev[i])
            if np.size(ind)>0:
                cfreq = len(ind[0])/len(ght)
            else:
                cfreq = 0.00001

            fftag=ftag+"ProbEvents_GEFShindcast2_"+repr(np.round(fqlev[i],2))
            true_binary = (gefs_hindcast_t > qlev[i]).astype(int)
            reliability_curve(true_binary,prob_results,cfreq,mlabels,nbins,pmax,fftag)

        except:
		    print("  - Did not compute Reliability Curve for "+repr(qlev[i]))
        else:
            print("  - Reliability Curve for "+repr(qlev[i])+" ok")

    plt.close('all')
    print(" Reliability Curve, OK")

