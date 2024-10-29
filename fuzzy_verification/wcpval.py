#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
wcpval.py

VERSION AND LAST UPDATE:
 v1.0  10/30/2023
 v1.1  11/03/2023
 v1.2  05/21/2024

PURPOSE:
 Validation of long-term probabilistic wave forecasts using 
 fuzzy verification. This can be used to assess and optimize
 the parameters of operational probability maps.
 The probability maps are based on GEFSv12 and have been validated against
  NDBC buoys, GDAS, and GEFS hindcasts by appending 24-hour segments 
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
 10/30/2023: Ricardo M. Campos, first version, named wprob.py
 11/03/2023: Ricardo M. Campos, fuzzy verification including categorical 
  and probabilistic assessments.
 05/21/2024: Ricardo M. Campos, renamed to wcpval.py, new reliability curve plot

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

import matplotlib
# matplotlib.use('Agg')
from matplotlib.dates import DateFormatter
import netCDF4 as nc
import numpy as np
import pandas as pd
import pickle
import sys
from sklearn.metrics import roc_curve, auc
from sklearn.calibration import calibration_curve
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import properscoring as ps
from datetime import datetime
import yaml
import warnings; warnings.filterwarnings("ignore")

sl=13
matplotlib.rcParams.update({'font.size': sl}); plt.rc('font', size=sl) 
matplotlib.rc('xtick', labelsize=sl); matplotlib.rc('ytick', labelsize=sl); matplotlib.rcParams.update({'font.size': sl})

# READ DATA
def read_data(wlist,bid,ltime1,ltime2,wvar):

    lstw=int(len(wlist))
    f=nc.Dataset(wlist[0])    
    stations = np.array(f.variables['buoyID'][bid]).astype('str')
    ensm = f.variables['ensemble_member'][:]
    latm = f.variables['lat'][0,:]; lonm = f.variables['lon'][0,:]
    indlat = int(np.floor(len(latm)/2)); indlon = int(np.floor(len(lonm)/2))
    auxt = np.array(f.variables['time'][:])
    # week 2 forecast
    axt = np.array(auxt-auxt.min())/3600.
    indt = np.where( (axt>=ltime1*24) & (axt<=ltime2*24) )[0]
    f.close(); del f, auxt, axt

    # read data
    for i in range(0,lstw):

        f=nc.Dataset(wlist[i])    

        if i==0:
            ctime = np.zeros((len(wlist)),'double')*np.nan
            ndbc = np.zeros((len(wlist),len(bid),len(indt)),'f')*np.nan
            fshape = f.variables[wvar+"_gefs_hindcast"].shape
            gefs_hindcast = np.zeros((len(wlist),len(bid),len(indt),fshape[2],fshape[3],fshape[4]),'f')*np.nan
            gefs_forecast = np.zeros((len(wlist),len(bid),len(indt),fshape[2],fshape[3],fshape[4]),'f')*np.nan
            del fshape

        ctime[i] = np.double(f.variables['time'][0])
        ndbc[i,:,:] = np.array(f.variables[wvar+"_ndbc"][bid,indt])
        gefs_hindcast[i,:,:,:,:,:] = np.array(f.variables[wvar+"_gefs_hindcast"][bid,indt,:,:,:])
        gefs_forecast[i,:,:,:,:,:] = np.array(f.variables[wvar+"_gefs_forecast"][bid,indt,:,:,:])

        f.close(); del f
        print(" read file "+repr(i)+" of "+repr(lstw-1))

    cdate = [datetime.utcfromtimestamp(ts) for ts in ctime]

    gefsdata = {'cdate': cdate,'ctime': ctime, 'stations': stations, 'lft': len(indt), 'ensm': ensm,
        'latm': latm, 'lonm': lonm, 'indlat': indlat, 'indlon': indlon, 
        'indlat': indlat, 'indlon': indlon, 'indt': indt,
        'ndbc': ndbc, 'gefs_hindcast': gefs_hindcast,'gefs_forecast': gefs_forecast}

    return gefsdata


# Function to check np.sort using multidimensional array
def tsort(model):

    fmodel=np.copy(model)
    if len(fmodel.shape) == 4:
        for i in range(0,model.shape[1]):
            for j in range(0,model.shape[2]):
                for k in range(0,model.shape[3]):
                    fmodel[:,i,j,k] = np.sort(model[:,i,j,k])

    elif len(fmodel.shape) == 3:
        for i in range(0,model.shape[1]):
            for j in range(0,model.shape[2]):
                fmodel[:,i,j] = np.sort(model[:,i,j])

    elif len(fmodel.shape) == 2:
        for i in range(0,model.shape[1]):
            fmodel[:,i] = np.sort(model[:,i])

    elif len(fmodel.shape) == 2:
        fmodel = np.sort(model[:])

    return fmodel


def nmaxsel(wdata,nmax):

    lstw = int(wdata.shape[0])
    wdata_tmax = np.mean(wdata,axis=1)
    for i in range(0,lstw):
        ind = np.where(wdata[i,:]>-999.)
        if np.size(ind)>int(np.floor(wdata.shape[1]/2)):
            ind = ind[0]
            wdata_tmax[i] = np.nanmean(np.sort(wdata[i,ind])[-nmax::])
        else:
            print(" nmaxsel(): forecast time series incomplete "+repr(i))
            wdata_tmax[i] = np.nan

        del ind

    return wdata_tmax


def probforecast(nmax,gspws,spctl,mod_forecast,qlev):

    lstw=int(mod_forecast.shape[0])
    nensm=int(mod_forecast.shape[2])
    indlat = int(np.floor(mod_forecast.shape[3]/2)); indlon = int(np.floor(mod_forecast.shape[4]/2))

    prob_forecast = np.zeros((lstw,len(qlev)),'f')*np.nan
    for i in range(0,lstw):
        fmod=np.array(mod_forecast[i,:,:,(indlat-gspws):(indlat+gspws+1),:][:,:,:,(indlon-gspws):(indlon+gspws+1)])
        fmod=tsort(fmod)
        fmod=np.copy(fmod[-nmax::,:,:,:])
        fmod=np.array(fmod.reshape(nmax*nensm,fmod.shape[2],fmod.shape[3]))
        fmod=np.array(fmod.reshape(nmax*nensm,fmod.shape[1]*fmod.shape[2]))
        fmod=np.sort(fmod,axis=1)
        ind=np.where(np.mean(fmod,axis=0)>=0.)
        if np.size(ind)>0:
            fmod=np.array(fmod[:,ind[0]])            
            fmod=np.array(fmod[:,int(np.floor(fmod.shape[1]*(spctl/100)))::])
            fmod=fmod.reshape(fmod.shape[0]*fmod.shape[1])
            for j in range(0,len(qlev)):
                prob_forecast[i,j] = np.size(fmod[fmod>qlev[j]]) / np.size(fmod[fmod>0.])

        if i==0:
            fmod_result = np.array([fmod])
        else:
            fmod_result = np.append(fmod_result,np.array([fmod]),axis=0)

    return prob_forecast, fmod_result


def prob_validation(nmax,spws,gspws,spctl,cdate,prob_u10_gefs_forecast,fmod_result_u10,prob_hs_gefs_forecast,fmod_result_hs,u10_obs_tmax,hs_obs_tmax,qlev_wnd,qlev_hs,plevels,ftag):

    # Binary Categorical array
    # U10
    fprob_u10_gefs_forecast = probforecast_binary(prob_u10_gefs_forecast,qlev_wnd,plevels)
    # Hs
    fprob_hs_gefs_forecast = probforecast_binary(prob_hs_gefs_forecast,qlev_hs,plevels)

    print(" Prob Forecast and Binary array, OK")

    lplev=int(len(plevels)-1)

    # ----------- Error metrics: TP, FP, TN, FN, CSI;  Brier; ROC_AUC ---------------
    ceval_gefs_u10 = np.zeros((len(qlev_wnd),lplev,3),'f')*np.nan
    for i in range(0,len(qlev_wnd)):
        for j in range(0,lplev):
            bresult = categorical_bin_eval(fprob_u10_gefs_forecast[:,i,j],u10_obs_tmax,qlev_wnd[i])
            ceval_gefs_u10[i,j,0] = float(bresult['POD'])
            ceval_gefs_u10[i,j,1] = float(bresult['FAR'])
            ceval_gefs_u10[i,j,2] = float(bresult['CSI'])
            del bresult

    ceval_gefs_hs = np.zeros((len(qlev_hs),lplev,3),'f')*np.nan
    for i in range(0,len(qlev_hs)):
        for j in range(0,lplev):
            bresult = categorical_bin_eval(fprob_hs_gefs_forecast[:,i,j],hs_obs_tmax,qlev_hs[i])
            ceval_gefs_hs[i,j,0] = float(bresult['POD'])
            ceval_gefs_hs[i,j,1] = float(bresult['FAR'])
            ceval_gefs_hs[i,j,2] = float(bresult['CSI'])
            del bresult

    # --- ROC Curve ---
    # U10
    froc_u10=[]
    for i in range(0,len(qlev_wnd)):
        fftag=ftag+"ProbEvents_U10_"+repr(np.round(qlev_wnd[i],2))
        true_binary = (u10_obs_tmax > qlev_wnd[i]).astype(int)
        roc = roc_plot(true_binary,prob_u10_gefs_forecast[:,i],fftag)
        froc_u10 = np.append(froc_u10,roc); del true_binary, fftag, roc

    # Hs
    froc_hs=[]
    for i in range(0,len(qlev_hs)):
        fftag=ftag+"ProbEvents_Hs_"+repr(np.round(qlev_hs[i],2))
        true_binary = (hs_obs_tmax > qlev_hs[i]).astype(int)
        roc = roc_plot(true_binary,prob_hs_gefs_forecast[:,i],fftag)
        froc_hs = np.append(froc_hs,roc); del true_binary, fftag, roc

    print(" Validation: Categorical OK")
    # --- Brier Score ---
    # U10
    fbriers_u10=[]
    for i in range(0,len(qlev_wnd)):    
        fftag=ftag+"ProbEvents_U10_"+repr(np.round(qlev_wnd[i],2))
        briers = brier_score(prob_u10_gefs_forecast[:,i],u10_obs_tmax,qlev_wnd[i],cdate,fftag)
        fbriers_u10 = np.append(fbriers_u10,briers)
        del briers, fftag

    # Hs
    fbriers_hs=[]
    for i in range(0,len(qlev_hs)):    
        fftag=ftag+"ProbEvents_Hs_"+repr(np.round(qlev_hs[i],2))
        briers = brier_score(prob_hs_gefs_forecast[:,i],hs_obs_tmax,qlev_hs[i],cdate,fftag)
        fbriers_hs = np.append(fbriers_hs,briers)
        del briers, fftag

    print(" Validation: Brier Score OK")

    # --- CRPS ---
    # U10 
    fftag=ftag+"ProbEvents_U10"
    crps_u10, mean_crps_u10 = crps(fmod_result_u10, u10_obs_tmax, cdate, fftag)
    del fftag
    # Hs
    fftag=ftag+"ProbEvents_hs"
    crps_hs, mean_crps_hs = crps(fmod_result_hs, hs_obs_tmax, cdate, fftag)
    del fftag

    plt.close('all')
    print(" Validation: CRPS OK")

    # ================================================

    # Save results
    stat_results = {
    'ceval_gefs_hs': ceval_gefs_hs,
    'ceval_gefs_u10': ceval_gefs_u10,
    'froc_hs': froc_hs,
    'froc_u10': froc_u10,
    'fbriers_u10': fbriers_u10,
    'fbriers_hs': fbriers_hs,
    'crps_u10': crps_u10,
    'mean_crps_u10': mean_crps_u10,
    'crps_hs': crps_hs,
    'mean_crps_hs': mean_crps_hs
    }

    # Save the dictionary to a pickle file
    with open(ftag+'STAT.RESULTS.pkl', 'wb') as f:
        pickle.dump(stat_results, f)

    return stat_results
    print(" Results saved in the pickle file STAT.RESULTS.pkl")

    del u10_obs_tmax, hs_obs_tmax, prob_u10_gefs_forecast, prob_hs_gefs_forecast, fprob_u10_gefs_forecast, fprob_hs_gefs_forecast
    del stat_results, ceval_gefs_u10, ceval_gefs_hs, froc_u10, froc_hs, fbriers_u10, fbriers_hs, crps_u10, mean_crps_u10, crps_hs, mean_crps_hs

    print(" - Done Optmz_nmax"+str(int(nmax))+"_spws"+str(int(spws*100)).zfill(3)+"_spctl"+str(int(spctl)))



# Binary Categorical array
def probforecast_binary(prob_forcast,qlev,plevels):

    lstw = int(prob_forcast.shape[0])
    lplev = int(len(plevels)-1)
    fprob_forecast = np.zeros((lstw,len(qlev),lplev),'f')*np.nan

    for i in range(0,lstw):
        for j in range(0,len(qlev)):
            for k in range(0,lplev):
                if prob_forcast[i,j] >= plevels[k]:
                    fprob_forecast[i,j,k] = int(1)
                else:
                    fprob_forecast[i,j,k] = int(0)

    return fprob_forecast


def categorical_bin_eval(forecasted,obs,threshold):

    ind = np.where((forecasted>-999.)&(obs>-999.))
    if np.size(ind)>1:
        forecasted=forecasted[ind[0]]
        obs=obs[ind[0]]

        SIZE = int(len(obs))
        EVENTS = int(np.size(np.where((obs >= threshold)==True)))

        # Calculate True Positives (TP), False Positives (FP), False Negatives (FN), and True Negatives (TN)
        if np.any((obs >= threshold)==True):
            TP = np.sum((forecasted == int(1)) & (obs >= threshold))
            FP = np.sum((forecasted == int(1)) & (obs < threshold))
            FN = np.sum((forecasted == int(0)) & (obs >= threshold))
            TN = np.sum((forecasted == int(0)) & (obs < threshold))

            # POD (Probability of Detection) also name as True Positive Rate (TPR) or Sensitivity
            if (TP + FN) == 0:
                POD = np.nan
            else:
                POD = np.round((TP / (TP + FN)),3)

            # FAR (False Alarm Ratio)
            if (FP + TP)==0:
                FAR = 0
            else:
                FAR = np.round((FP / (FP + TP)),3)

            # Critical Success Index (CSI)
            if (TP + FN + FP)==0:
                CSI = np.nan
            else:
                CSI = np.round((TP / (TP + FN + FP)),3)

            # ETS (Equitable Threat Score)
            ETS = np.round(((TP - (TP + FN) * (TP + FP) / (TP + FN + FP + TN)) / (TP + FN + FP - (TP + FN) * (TP + FP) / (TP + FN + FP + TN))) ,3)
        else:
            TP=np.nan; FP=np.nan; FN=np.nan; TN=np.nan
            POD=np.nan; FAR=np.nan; CSI=np.nan; ETS=np.nan

        result = {'SIZE':SIZE, 'EVENTS':EVENTS,
            'TP':TP, 'FP':FP, 'FN':FN, 'TN':TN,
            'POD':POD, 'FAR':FAR, 'CSI':CSI, 'ETS': ETS}

        return result

def brier_score(prob_forecast,obs,threshold,cdate,ftag):

    # Brier Score
    true_binary = (obs > threshold).astype(int)
    briers=0
    for i in range(0,len(true_binary)):
        briers = briers+(true_binary[i]-prob_forecast[i])**2

    briers = briers/len(true_binary)

    # Plot Prob and Brier Score Scheme
    fig1, ax = plt.subplots(figsize=(7, 4))
    ax.plot(cdate, true_binary, color='dimgray', marker='o', linestyle='', linewidth=2., label='Obs(Event)', zorder=2)
    # ax.plot(cdate, true_binary, color='black', marker='.', linestyle='', linewidth=1., zorder=2)
    ax.plot(cdate, prob_forecast, color='firebrick', marker='*', linestyle='', linewidth=2., label='ProbModel', zorder=2)
    ax.set_xlim(cdate[0], cdate[-1]); ax.set_ylim(-0.04, 1.04)
    ax.xaxis.set_major_formatter(DateFormatter('%b-%Y'))
    ax.fmt_xdata = DateFormatter('%b-%Y')
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(sl - 6)

    ax.set_xlabel("Time")
    ax.set_ylabel("Probability")
    ax.legend(loc='center left',fontsize=sl-5)
    plt.grid(c='grey', ls='--', alpha=0.3, zorder=1)
    plt.savefig(ftag+"_BrierScore.png", dpi=200, facecolor='w', edgecolor='w', orientation='portrait',
        format='png', bbox_inches='tight', pad_inches=0.1)

    plt.close(fig1)
    return float(briers)


def crps(fmod,obs,cdate,ftag):

    # CRPS
    crps = ps.crps_ensemble(obs,fmod)
    fcrps_mean = np.nanmean(crps)

    # Plot
    fig1, ax = plt.subplots(figsize=(7, 4))
    ax.vlines(cdate, 0, obs/np.nanmax(obs), colors='silver', linewidth=1., zorder=1)
    ax.plot(cdate, obs/np.nanmax(obs), color='grey', marker='.', linestyle='', linewidth=2., label='Obs (normalized)', zorder=2)
    ax.plot(cdate, crps, color='blue', marker='.', linestyle='', linewidth=2., label='CRPS', zorder=3)
    ax.set_xlim(cdate[0], cdate[-1]); ax.set_ylim(-0.04, np.nanmax(crps)+0.1)
    ax.xaxis.set_major_formatter(DateFormatter('%b-%Y'))
    ax.fmt_xdata = DateFormatter('%b-%Y')
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(sl - 6)

    ax.set_xlabel("Time")
    ax.set_ylabel("CRPS")
    ax.legend(loc='center left',fontsize=sl-5)
    plt.grid(c='grey', ls='--', alpha=0.3, zorder=1)
    plt.savefig(ftag+"_CRPS.png", dpi=200, facecolor='w', edgecolor='w', orientation='portrait',
        format='png', bbox_inches='tight', pad_inches=0.1)

    plt.close(fig1)
    return np.array(crps).astype('float'), float(fcrps_mean)
    del crps, fcrps_mean


def roc_plot(true_binary,predicted_probabilities,ftag):

    # Calculate ROC curve
    fpr, tpr, _ = roc_curve(true_binary, predicted_probabilities)
    # Calculate AUC (Area Under the ROC Curve)
    roc_auc = auc(fpr, tpr)

    if len(ftag)>0:
        # Plot the ROC curve
        fig1, ax = plt.subplots(figsize=(5,4.5))
        ax.plot([0, 1], [0, 1], color='dimgray', lw=2, linestyle='--', alpha=0.8, zorder=1)
        ax.plot(fpr, tpr, color='firebrick', lw=1, marker='.', alpha=0.8, zorder=2)
        ax.plot(fpr, tpr, color='firebrick', lw=2, label=f'ROC curve (AUC = {roc_auc:.2f})', zorder=3)
        ax.set_xlim(-0.01, 1.01); ax.set_ylim(-0.01, 1.01) 
        ax.set_xlabel('False Positive Rate')
        ax.set_ylabel('True Positive Rate')
        ax.legend(loc='lower right',fontsize=sl-3)
        plt.title('ROC Curve')
        plt.grid(c='grey', ls='--', alpha=0.3, zorder=1)
        plt.savefig(ftag+"_ROC.png", dpi=200, facecolor='w', edgecolor='w', orientation='portrait',
            format='png', bbox_inches='tight', pad_inches=0.1)

        plt.close(fig1)

    return float(roc_auc)

def reliability_curve(true_binary,predicted_probabilities,ftag):

    #ind = np.where( (prob_forecast>-999) & (obs>-999) )[0]
    #reldiag = verification.reldiag_init(threshold)
    #verification.reldiag_accum(reldiag, prob_forecast[ind],obs[ind])
    #fig, ax = plt.subplots()
    #verification.plot_reldiag(reldiag, ax)

    ry, rx = calibration_curve(true_binary, predicted_probabilities, n_bins=10)

    fig1, ax = plt.subplots(figsize=(5,4.5))
    ax.plot(rx, ry, marker='s', linewidth=1, color='navy', zorder=3)
    line = mlines.Line2D([0, 1], [0, 1], color='dimgrey',linestyle='--', zorder=2)
    transform = ax.transAxes; line.set_transform(transform); ax.add_line(line)
    ax.set_xlabel('Forecast probability')
    ax.set_ylabel('Observed relative frequency')
    plt.grid(c='grey', ls='--', alpha=0.3, zorder=1)
    plt.ylim(ymax = 1.01, ymin = -0.01)
    plt.xlim(xmax = 1.01, xmin = -0.01)
    plt.tight_layout()
    # ax.set_title("Reliability diagram" )
    plt.savefig(ftag+"_ReliabilityDiagram.png", dpi=200, facecolor='w', edgecolor='w', orientation='portrait',
        format='png', bbox_inches='tight', pad_inches=0.1)

    plt.close(fig1)



if __name__ == "__main__":

    # select one point
    # station = '46006'
    station=str(sys.argv[1])
    # Forecast Lead Time (Day) and intervall
    # ltime1=7
    ltime1=int(sys.argv[2])
    # ltime2=14
    ltime2=int(sys.argv[3])
    # opath="/home/ricardo/cimas/analysis/4postproc/output"
    opath=str(sys.argv[4])

    # file tag to same output fils
    ftag=opath+"/Validation_"+station+"_"

    # ---- Read statistical parameters (exaclty the same as the operational config file) -----
    print(" "); print(" Reading yaml configuration file ...")
    with open('probmaps_gefs.yaml', 'r') as file:
        wconfig = yaml.safe_load(file)

    qlev_hs = np.array(wconfig['qlev_hs']).astype('float'); fqlev_hs=qlev_hs
    fqlev_wnd = np.array(wconfig['qlev_wnd']).astype('float')
    qlev_wnd = fqlev_wnd/1.94 # m/s
    plevels = np.array(wconfig['plevels']).astype('float')
    nmax = int(wconfig['nmax'])

    # -------------------------------------

    # READ DATA
    # list of netcdf files generated with buildfuzzydataset.py (GEFS, GDAS, and NDBC buoy)
    # ls -d $PWD/*.nc > list.txt &
    wlist=np.atleast_1d(np.loadtxt('list.txt',dtype=str)) 
    lstw=int(len(wlist)); lplev=int(len(plevels)-1)

    cdate,ctime,u10_ndbc,hs_ndbc,u10_gdas,hs_gdas,u10_gefs_hindcast,hs_gefs_forecast = read_data(wlist,station)

    # Preserve original data and Reshape
    ru10_ndbc = np.array(u10_ndbc.reshape(lstw*len(indt)))
    rhs_ndbc = np.array(hs_ndbc.reshape(lstw*len(indt)))
    ru10_gdas = np.array(u10_gdas.reshape((lstw*len(indt),len(latm),len(lonm))))
    rhs_gdas = np.array(hs_gdas.reshape((lstw*len(indt),len(latm),len(lonm))))
    ru10_gefs_hindcast = np.array(u10_gefs_hindcast.reshape((lstw*len(indt),len(ensm),len(latm),len(lonm))))
    rhs_gefs_hindcast = np.array(hs_gefs_hindcast.reshape((lstw*len(indt),len(ensm),len(latm),len(lonm))))
    ru10_gefs_forecast = np.array(u10_gefs_forecast.reshape((lstw*len(indt),len(ensm),len(latm),len(lonm))))
    rhs_gefs_forecast = np.array(hs_gefs_forecast.reshape((lstw*len(indt),len(ensm),len(latm),len(lonm))))

    # ===== Mimic the operational probability maps methodology ===== 

    u10_ndbc_t = nmaxsel(u10_ndbc,nmax)
    hs_ndbc_t = nmaxsel(hs_ndbc,nmax)

    #  -- Probabilistic Forecast array --
    
    nmax = int(wconfig['nmax']); spws = int(wconfig['spws']); spctl = float(wconfig['spctl'])
    gspws=int(np.floor(spws/np.diff(latm).mean())/2)
    # U10 
    prob_u10_gefs_forecast, fmod_result_u10 = probforecast(nmax,gspws,spctl,u10_gefs_forecast,qlev_wnd)
    # Hs
    prob_hs_gefs_forecast, fmod_result_hs = probforecast(nmax,gspws,spctl,hs_gefs_forecast,qlev_hs)


    #  ++++++++ VALIDATION +++++++++++++++

    # ===== Categorical approach, verification of binary events =====

    # Binary Categorical array
    # U10
    fprob_u10_gefs_forecast = probforecast_binary(prob_u10_gefs_forecast,qlev_wnd,plevels)
    # Hs
    fprob_hs_gefs_forecast = probforecast_binary(prob_hs_gefs_forecast,qlev_hs,plevels)

    # Validation scores
    # U10
    for i in range(0,len(qlev_wnd)):
        for j in range(0,lplev):
            ceval_gefs_ndbc = categorical_bin_eval(fprob_u10_gefs_forecast[:,i,j],u10_ndbc_t,qlev_wnd[i])
            ceval_gefs_gdas = categorical_bin_eval(fprob_u10_gefs_forecast[:,i,j],u10_gdas_t,qlev_wnd[i])
            ceval_gefs_ghnd = categorical_bin_eval(fprob_u10_gefs_forecast[:,i,j],u10_gefs_hindcast_t,qlev_wnd[i])
            dict_list = [
                {'Name': 'gefs_ndbc', **ceval_gefs_ndbc},
                {'Name': 'gefs_gdas', **ceval_gefs_gdas},
                {'Name': 'gefs_gefsHindcast', **ceval_gefs_ghnd}
            ]

            # Convert the list of dictionaries into a pandas DataFrame
            df = pd.DataFrame(dict_list)
            df.to_csv(ftag+"CategoricalVal_U10_Lev"+repr(np.round(fqlev_wnd[i],2))+".csv",sep='\t', index=False)

    # Hs
    for i in range(0,len(qlev_hs)):
        for j in range(0,lplev):
            ceval_gefs_ndbc = categorical_bin_eval(fprob_hs_gefs_forecast[:,i,j],hs_ndbc_t,qlev_hs[i])
            ceval_gefs_gdas = categorical_bin_eval(fprob_hs_gefs_forecast[:,i,j],hs_gdas_t,qlev_hs[i])
            ceval_gefs_ghnd = categorical_bin_eval(fprob_hs_gefs_forecast[:,i,j],hs_gefs_hindcast_t,qlev_hs[i])
            dict_list = [
                {'Name': 'gefs_ndbc', **ceval_gefs_ndbc},
                {'Name': 'gefs_gdas', **ceval_gefs_gdas},
                {'Name': 'gefs_gefsHindcast', **ceval_gefs_ghnd}
            ]

            # Convert the list of dictionaries into a pandas DataFrame
            df = pd.DataFrame(dict_list)
            df.to_csv(ftag+"CategoricalVal_Hs_Lev"+repr(np.round(fqlev_hs[i],2))+"_Prob"+repr(np.round(plevels[j],2))+".csv",sep='\t', index=False)


    # Probabilistic approach, verification of event frequency
    # BS, CRPS, ROC, reliability curve

    onames = ['gefs_ndbc', 'gefs_gdas', 'gefs_gefsHindcast']

    plt.close('all')
    # --- Brier Score ---
    # U10
    for i in range(0,len(qlev_wnd)):
        fbriers=[]
        fftag=ftag+"ProbEvents_NDBC_U10_"+repr(np.round(fqlev_wnd[i],2))
        briers = brier_score(prob_u10_gefs_forecast[:,i],u10_ndbc_t,qlev_wnd[i],cdate,fftag)
        fbriers = np.append(fbriers,briers)

        fftag=ftag+"ProbEvents_GDAS_U10_"+repr(np.round(fqlev_wnd[i],2))
        briers = brier_score(prob_u10_gefs_forecast[:,i],u10_gdas_t,qlev_wnd[i],cdate,fftag)
        fbriers = np.append(fbriers,briers)

        fftag=ftag+"ProbEvents_GEFShindcast_U10_"+repr(np.round(fqlev_wnd[i],2))
        briers = brier_score(prob_u10_gefs_forecast[:,i],u10_gefs_hindcast_t,qlev_wnd[i],cdate,fftag)
        fbriers = np.append(fbriers,briers)

        bdata = {'Name': onames, 'Result': np.round(fbriers,3)}
        df = pd.DataFrame(bdata)
        df.to_csv(ftag+"ProbabilisticVal_BrierScore_U10_Lev"+repr(np.round(fqlev_wnd[i],2))+".csv",sep='\t', index=False)

    # Hs
    for i in range(0,len(qlev_hs)):
        fbriers=[]
        fftag=ftag+"ProbEvents_NDBC_Hs_"+repr(np.round(fqlev_hs[i],2))
        briers = brier_score(prob_hs_gefs_forecast[:,i],hs_ndbc_t,qlev_hs[i],cdate,fftag)
        fbriers = np.append(fbriers,briers)

        fftag=ftag+"ProbEvents_GDAS_Hs_"+repr(np.round(fqlev_hs[i],2))
        briers = brier_score(prob_hs_gefs_forecast[:,i],hs_gdas_t,qlev_hs[i],cdate,fftag)
        fbriers = np.append(fbriers,briers)

        fftag=ftag+"ProbEvents_GEFShindcast_Hs_"+repr(np.round(fqlev_hs[i],2))
        briers = brier_score(prob_hs_gefs_forecast[:,i],hs_gefs_hindcast_t,qlev_hs[i],cdate,fftag)
        fbriers = np.append(fbriers,briers)

        bdata = {'Name': onames, 'Result': np.round(fbriers,3)}
        df = pd.DataFrame(bdata)
        df.to_csv(ftag+"ProbabilisticVal_BrierScore_Hs_Lev"+repr(np.round(fqlev_hs[i],2))+".csv",sep='\t', index=False)

    plt.close('all')
    # --- CRPS ---
    # U10 
    fcrps=[]
    fftag=ftag+"ProbEvents_NDBC_U10"
    briers, crps = brier_crps(fmod_result_u10,prob_u10_gefs_forecast[:,i],u10_ndbc_t,qlev_wnd[0],cdate,fftag)
    fcrps = np.append(fcrps,crps)

    fftag=ftag+"ProbEvents_GDAS_U10"
    briers, crps = brier_crps(fmod_result_u10,prob_u10_gefs_forecast[:,i],u10_gdas_t,qlev_wnd[0],cdate,fftag)
    fcrps = np.append(fcrps,crps)

    fftag=ftag+"ProbEvents_GEFShindcast_U10"
    briers, crps = brier_crps(fmod_result_u10,prob_u10_gefs_forecast[:,i],u10_gefs_hindcast_t,qlev_wnd[0],cdate,fftag)
    fcrps = np.append(fcrps,crps)

    bdata = {'Name': onames, 'Result': np.round(fcrps,3)}
    df = pd.DataFrame(bdata)
    df.to_csv(ftag+"ProbabilisticVal_CRPS_U10.csv",sep='\t', index=False)

    # Hs 
    fcrps=[]
    fftag=ftag+"ProbEvents_NDBC_Hs"
    briers, crps = brier_crps(fmod_result_hs,prob_hs_gefs_forecast[:,i],hs_ndbc_t,qlev_hs[0],cdate,fftag)
    fcrps = np.append(fcrps,crps)

    fftag=ftag+"ProbEvents_GDAS_Hs"
    briers, crps = brier_crps(fmod_result_hs,prob_hs_gefs_forecast[:,i],hs_gdas_t,qlev_hs[0],cdate,fftag)
    fcrps = np.append(fcrps,crps)

    fftag=ftag+"ProbEvents_GEFShindcast_Hs"
    briers, crps = brier_crps(fmod_result_hs,prob_hs_gefs_forecast[:,i],hs_gefs_hindcast_t,qlev_hs[0],cdate,fftag)
    fcrps = np.append(fcrps,crps)

    bdata = {'Name': onames, 'Result': np.round(fcrps,3)}
    df = pd.DataFrame(bdata)
    df.to_csv(ftag+"ProbabilisticVal_CRPS_Hs.csv",sep='\t', index=False)

    plt.close('all')
    # --- ROC Curve ---
    # U10
    for i in range(0,len(qlev_wnd)):
        froc=[]
        fftag=ftag+"ProbEvents_NDBC_U10_"+repr(np.round(fqlev_wnd[i],2))
        true_binary = (u10_ndbc_t > qlev_wnd[i]).astype(int)
        roc = roc_plot(true_binary,prob_u10_gefs_forecast[:,i],fftag)
        froc = np.append(froc,roc); del true_binary, fftag

        fftag=ftag+"ProbEvents_GDAS_U10_"+repr(np.round(fqlev_wnd[i],2))
        true_binary = (u10_gdas_t > qlev_wnd[i]).astype(int)
        roc = roc_plot(true_binary,prob_u10_gefs_forecast[:,i],fftag)
        froc = np.append(froc,roc); del true_binary, fftag

        fftag=ftag+"ProbEvents_GEFShindcast_U10_"+repr(np.round(fqlev_wnd[i],2))
        true_binary = (u10_gefs_hindcast_t > qlev_wnd[i]).astype(int)
        roc = roc_plot(true_binary,prob_u10_gefs_forecast[:,i],fftag)
        froc = np.append(froc,roc); del true_binary, fftag

        bdata = {'Name': onames, 'Result': np.round(froc,3)}
        df = pd.DataFrame(bdata)
        df.to_csv(ftag+"ProbabilisticVal_ROCauc_U10_Lev"+repr(np.round(fqlev_wnd[i],2))+".csv",sep='\t', index=False)

    # Hs
    for i in range(0,len(qlev_hs)):
        froc=[]
        fftag=ftag+"ProbEvents_NDBC_Hs_"+repr(np.round(fqlev_hs[i],2))
        true_binary = (hs_ndbc_t > qlev_hs[i]).astype(int)
        roc = roc_plot(true_binary,prob_hs_gefs_forecast[:,i],fftag)
        froc = np.append(froc,roc); del true_binary, fftag

        fftag=ftag+"ProbEvents_GDAS_Hs_"+repr(np.round(fqlev_hs[i],2))
        true_binary = (hs_gdas_t > qlev_hs[i]).astype(int)
        roc = roc_plot(true_binary,prob_hs_gefs_forecast[:,i],fftag)
        froc = np.append(froc,roc); del true_binary, fftag

        fftag=ftag+"ProbEvents_GEFShindcast_Hs_"+repr(np.round(fqlev_hs[i],2))
        true_binary = (hs_gefs_hindcast_t > qlev_hs[i]).astype(int)
        roc = roc_plot(true_binary,prob_hs_gefs_forecast[:,i],fftag)
        froc = np.append(froc,roc); del true_binary, fftag

        bdata = {'Name': onames, 'Result': np.round(froc,3)}
        df = pd.DataFrame(bdata)
        df.to_csv(ftag+"ProbabilisticVal_ROCauc_Hs_Lev"+repr(np.round(fqlev_hs[i],2))+".csv",sep='\t', index=False)

    plt.close('all')
    # --- Reliability Curve ---
    # U10
#    for i in range(0,len(qlev_wnd)):
#        fftag=ftag+"ProbEvents_NDBC_U10_"+repr(np.round(fqlev_wnd[i],2))
#        reliability_curve(prob_u10_gefs_forecast[:,i],u10_ndbc_t,qlev_wnd[i],fftag)

#        fftag=ftag+"ProbEvents_GDAS_U10_"+repr(np.round(fqlev_wnd[i],2))
#        reliability_curve(prob_u10_gefs_forecast[:,i],u10_gdas_t,qlev_wnd[i],fftag)

#        fftag=ftag+"ProbEvents_GEFShindcast_U10_"+repr(np.round(fqlev_wnd[i],2))
#        reliability_curve(prob_u10_gefs_forecast[:,i],u10_gefs_hindcast_t,qlev_wnd[i],fftag)

    # Hs
#    for i in range(0,len(qlev_hs)):
#        fftag=ftag+"ProbEvents_NDBC_Hs_"+repr(np.round(fqlev_hs[i],2))
#        reliability_curve(prob_hs_gefs_forecast[:,i],hs_ndbc_t,qlev_hs[i],fftag)

#        fftag=ftag+"ProbEvents_GDAS_Hs_"+repr(np.round(fqlev_hs[i],2))
#        reliability_curve(prob_hs_gefs_forecast[:,i],hs_gdas_t,qlev_hs[i],fftag)

#        fftag=ftag+"ProbEvents_GEFShindcast_Hs_"+repr(np.round(fqlev_hs[i],2))
#        reliability_curve(prob_hs_gefs_forecast[:,i],hs_gefs_hindcast_t,qlev_hs[i],fftag)

    plt.close('all')

