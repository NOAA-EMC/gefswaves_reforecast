#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
fuzzy_verification_ProbMaps.py

VERSION AND LAST UPDATE:
 v1.0  10/30/2023
 v1.1  11/03/2023

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
 10/30/2023: Ricardo M. Campos, first version.
 11/03/2023: Ricardo M. Campos, fuzzy verification including categorical 
  and probabilistic assessments.

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

import matplotlib
matplotlib.use('Agg')
from matplotlib.dates import DateFormatter
import netCDF4 as nc
import numpy as np
import sys
import pandas as pd
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
import properscoring as ps
from datetime import datetime
from pysteps import verification
import yaml
from pvalstats import ModelObsPlot
import warnings; warnings.filterwarnings("ignore")

sl=13
matplotlib.rcParams.update({'font.size': sl}); plt.rc('font', size=sl) 
matplotlib.rc('xtick', labelsize=sl); matplotlib.rc('ytick', labelsize=sl); matplotlib.rcParams.update({'font.size': sl})

# READ DATA
def read_data(wlist,station,ltime1,ltime2,lstw):

    f=nc.Dataset(wlist[0])    
    st = np.where(f.variables['buoyID'][:]==station)[0][0]
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
            ctime=np.array([f.variables['time'][0]]).astype('double')
            u10_ndbc = np.array([f.variables['u10_ndbc'][st,indt]])
            hs_ndbc = np.array([f.variables['hs_ndbc'][st,indt]])
            u10_gdas = np.array([f.variables['u10_gdas'][st,indt,:,:]])
            hs_gdas = np.array([f.variables['hs_gdas'][st,indt,:,:]])
            u10_gefs_hindcast = np.array([f.variables['u10_gefs_hindcast'][st,indt,:,:,:]])
            hs_gefs_hindcast = np.array([f.variables['hs_gefs_hindcast'][st,indt,:,:,:]])
            u10_gefs_forecast = np.array([f.variables['u10_gefs_forecast'][st,indt,:,:,:]])
            hs_gefs_forecast = np.array([f.variables['hs_gefs_forecast'][st,indt,:,:,:]])
        else:
            ctime = np.append(ctime, np.array([f.variables['time'][0]]).astype('double'))
            u10_ndbc = np.append(u10_ndbc,np.array([f.variables['u10_ndbc'][st,indt]]),axis=0)
            hs_ndbc = np.append(hs_ndbc,np.array([f.variables['hs_ndbc'][st,indt]]),axis=0)
            u10_gdas = np.append(u10_gdas,np.array([f.variables['u10_gdas'][st,indt,:,:]]),axis=0)
            hs_gdas = np.append(hs_gdas,np.array([f.variables['hs_gdas'][st,indt,:,:]]),axis=0)
            u10_gefs_hindcast = np.append(u10_gefs_hindcast,np.array([f.variables['u10_gefs_hindcast'][st,indt,:,:,:]]),axis=0)
            hs_gefs_hindcast = np.append(hs_gefs_hindcast,np.array([f.variables['hs_gefs_hindcast'][st,indt,:,:,:]]),axis=0)
            u10_gefs_forecast = np.append(u10_gefs_forecast,np.array([f.variables['u10_gefs_forecast'][st,indt,:,:,:]]),axis=0)
            hs_gefs_forecast = np.append(hs_gefs_forecast,np.array([f.variables['hs_gefs_forecast'][st,indt,:,:,:]]),axis=0)

        f.close(); del f
        print(" read file "+repr(i)+" of "+repr(lstw-1))

    cdate = [datetime.utcfromtimestamp(ts) for ts in ctime]

    return cdate,ctime,st,ensm,latm,lonm,indlat,indlon,u10_ndbc,hs_ndbc,u10_gdas,hs_gdas,u10_gefs_hindcast,hs_gefs_hindcast,u10_gefs_forecast,hs_gefs_forecast,indt

def categorical_eval(forecasted,obs,threshold):

    ind = np.where((forecasted>-999.)&(obs>-999.))
    if np.size(ind)>1:
        forecasted=forecasted[ind[0]]
        obs=obs[ind[0]]

        SIZE = int(len(obs))
        EVENTS = int(np.size(np.where((obs >= threshold)==True)))

        # Calculate True Positives (TP), False Positives (FP), False Negatives (FN), and True Negatives (TN)
        TP = np.sum((forecasted >= threshold) & (obs >= threshold))
        FP = np.sum((forecasted >= threshold) & (obs < threshold))
        FN = np.sum((forecasted < threshold) & (obs >= threshold))
        TN = np.sum((forecasted < threshold) & (obs < threshold))

        # POD (Probability of Detection) also name as True Positive Rate (TPR) or Sensitivity
        POD = np.round((TP / (TP + FN)),3)
        # FAR (False Alarm Ratio)
        FAR = np.round((FP / (FP + TP)),3)
        # Critical Success Index (CSI)
        CSI = np.round((TP / (TP + FN + FP)),3)
        # ETS (Equitable Threat Score)
        ETS = np.round(((TP - (TP + FN) * (TP + FP) / (TP + FN + FP + TN)) / (TP + FN + FP - (TP + FN) * (TP + FP) / (TP + FN + FP + TN))) ,3)

        result = {'SIZE':SIZE, 'EVENTS':EVENTS,
            'TP':TP, 'FP':FP, 'FN':FN, 'TN':TN,
            'POD':POD, 'FAR':FAR, 'CSI':CSI, 'ETS': ETS}

        return result

def categorical_bin_eval(forecasted,obs,threshold):

    ind = np.where((forecasted>-999.)&(obs>-999.))
    if np.size(ind)>1:
        forecasted=forecasted[ind[0]]
        obs=obs[ind[0]]

        SIZE = int(len(obs))
        EVENTS = int(np.size(np.where((obs >= threshold)==True)))

        # Calculate True Positives (TP), False Positives (FP), False Negatives (FN), and True Negatives (TN)
        TP = np.sum((forecasted == int(1)) & (obs >= threshold))
        FP = np.sum((forecasted == int(1)) & (obs < threshold))
        FN = np.sum((forecasted == int(0)) & (obs >= threshold))
        TN = np.sum((forecasted == int(0)) & (obs < threshold))

        # POD (Probability of Detection) also name as True Positive Rate (TPR) or Sensitivity
        POD = np.round((TP / (TP + FN)),3)
        # FAR (False Alarm Ratio)
        FAR = np.round((FP / (FP + TP)),3)
        # Critical Success Index (CSI)
        CSI = np.round((TP / (TP + FN + FP)),3)
        # ETS (Equitable Threat Score)
        ETS = np.round(((TP - (TP + FN) * (TP + FP) / (TP + FN + FP + TN)) / (TP + FN + FP - (TP + FN) * (TP + FP) / (TP + FN + FP + TN))) ,3)

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
    fig1, ax = plt.subplots(figsize=(6, 4))
    ax.plot(cdate, true_binary, color='dimgray', marker='o', linestyle='', linewidth=2., label='Obs(Event)', zorder=2)
    # ax.plot(cdate, true_binary, color='black', marker='.', linestyle='', linewidth=1., zorder=2)
    ax.plot(cdate, prob_forecast, color='firebrick', marker='*', linestyle='', linewidth=2., label='ProbModel', zorder=2)
    ax.set_xlim(cdate[0], cdate[-1]); ax.set_ylim(-0.04, 1.04)
    ax.xaxis.set_major_formatter(DateFormatter('%b%d'))
    ax.fmt_xdata = DateFormatter('%b%d')
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(sl - 4)

    ax.set_xlabel("Time")
    ax.set_ylabel("Probability")
    ax.legend(loc='center left',fontsize=sl-5)
    plt.grid(c='grey', ls='--', alpha=0.3, zorder=1)
    plt.savefig(ftag+"_BrierScore.png", dpi=200, facecolor='w', edgecolor='w', orientation='portrait',
        format='png', bbox_inches='tight', pad_inches=0.1)

    plt.close(fig1)
    return float(briers)


def brier_crps(fmod,prob_forecast,obs,threshold,cdate,ftag):

    # Brier Score
    true_binary = (obs > threshold).astype(int)
    briers=0
    for i in range(0,len(true_binary)):
        briers = briers+(true_binary[i]-prob_forecast[i])**2

    briers = briers/len(true_binary)

    # CRPS
    crps = ps.crps_ensemble(obs,fmod)
    fcrps = np.nanmean(crps)

    # Plot CRPS and Prob
    fig1, ax = plt.subplots(figsize=(6, 4))
    ax.plot(cdate, true_binary, color='dimgray', marker='o', linestyle='', linewidth=2., label='Obs(Event)', alpha=0.4, zorder=1)
    ax.plot(cdate, prob_forecast, color='firebrick', marker='*', linestyle='', linewidth=2., label='ProbModel', alpha=0.4, zorder=2)
    ax.plot(cdate, crps, color='navy', marker='.', linestyle='', linewidth=2., label='CRPS', zorder=3)
    ax.set_xlim(cdate[0], cdate[-1]); ax.set_ylim(-0.04, np.nanmax(crps)+0.1)
    ax.xaxis.set_major_formatter(DateFormatter('%b%d'))
    ax.fmt_xdata = DateFormatter('%b%d')
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(sl - 4)

    ax.set_xlabel("Time")
    ax.set_ylabel("Probability & CRPS")
    ax.legend(loc='center left',fontsize=sl-5)
    plt.grid(c='grey', ls='--', alpha=0.3, zorder=1)
    plt.savefig(ftag+"_CRPS.png", dpi=200, facecolor='w', edgecolor='w', orientation='portrait',
        format='png', bbox_inches='tight', pad_inches=0.1)

    plt.close(fig1)

    return float(briers),float(fcrps)

def roc_plot(true_binary,predicted_probabilities,ftag):

    # Calculate ROC curve
    fpr, tpr, _ = roc_curve(true_binary, predicted_probabilities)
    # Calculate AUC (Area Under the ROC Curve)
    roc_auc = auc(fpr, tpr)

    # Plot the ROC curve
    fig1, ax = plt.subplots(figsize=(6, 5))
    ax.plot([0, 1], [0, 1], color='dimgray', lw=2, linestyle='--', alpha=0.8, zorder=1)
    ax.plot(fpr, tpr, color='firebrick', lw=1, marker='.', alpha=0.8, zorder=2)
    ax.plot(fpr, tpr, color='firebrick', lw=2, label=f'ROC curve (AUC = {roc_auc:.2f})', zorder=3)
    ax.set_xlim(-0.01, 1.01); ax.set_ylim(-0.01, 1.01) 
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.legend(loc='lower right',fontsize=sl-5)
    plt.title('ROC Curve')
    plt.grid(c='grey', ls='--', alpha=0.3, zorder=1)
    plt.savefig(ftag+"_ROC.png", dpi=200, facecolor='w', edgecolor='w', orientation='portrait',
        format='png', bbox_inches='tight', pad_inches=0.1)

    plt.close(fig1)

    return float(roc_auc)

def reliability_curve(prob_forecast,obs,threshold,ftag):

    reldiag = verification.reldiag_init(threshold)
    verification.reldiag_accum(reldiag, prob_forecast,obs)
    fig, ax = plt.subplots()
    verification.plot_reldiag(reldiag, ax)
    ax.set_title("Reliability diagram" )
    plt.savefig(ftag+"_ReliabilityDiagram.png", dpi=200, facecolor='w', edgecolor='w', orientation='portrait',
        format='png', bbox_inches='tight', pad_inches=0.1)

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

    # -- Selec nmax --
    u10_ndbc_t = np.mean(u10_ndbc,axis=1)
    for i in range(0,lstw):
        ind = np.where(u10_ndbc[i,:]>-999.)
        if np.size(ind)>int(np.floor(len(indt)/2)) :
            ind = ind[0]
            u10_ndbc_t[i] = np.nanmean(np.sort(u10_ndbc[i,ind])[-nmax::])

        del ind

    hs_ndbc_t = np.mean(hs_ndbc,axis=1)
    for i in range(0,lstw):
        ind = np.where(hs_ndbc[i,:]>-999.)
        if np.size(ind)>int(np.floor(len(indt)/2)) :
            ind = ind[0]
            hs_ndbc_t[i] = np.nanmean(np.sort(hs_ndbc[i,ind])[-nmax::])

        del ind

    u10_gdas_t = np.mean(u10_gdas[:,:,indlat,indlon],axis=1)
    for i in range(0,lstw):
        ind = np.where(u10_gdas[i,:,indlat,indlon]>-999.)
        if np.size(ind)>int(np.floor(len(indt)/2)) :
            ind = ind[0]
            u10_gdas_t[i] = np.nanmean(np.sort(u10_gdas[i,ind,indlat,indlon])[-nmax::])

        del ind

    hs_gdas_t = np.mean(hs_gdas[:,:,indlat,indlon],axis=1)
    for i in range(0,lstw):
        ind = np.where(hs_gdas[i,:,indlat,indlon]>-999.)
        if np.size(ind)>int(np.floor(len(indt)/2)) :
            ind = ind[0]
            hs_gdas_t[i] = np.nanmean(np.sort(hs_gdas[i,ind,indlat,indlon])[-nmax::])

        del ind

    u10_gefs_hindcast_t = np.mean(u10_gefs_hindcast[:,:,:,indlat,indlon],axis=(1,2))
    for i in range(0,lstw):
        aux = np.nanmean(u10_gefs_hindcast[i,:,:,indlat,indlon],axis=1) # average ensemble members
        ind = np.where(aux>-999.)
        if np.size(ind)>int(np.floor(len(indt)/2)) :
            ind = ind[0]
            u10_gefs_hindcast_t[i] = np.nanmean(np.sort(aux[ind])[-nmax::])

        del ind,aux

    hs_gefs_hindcast_t = np.mean(hs_gefs_hindcast[:,:,:,indlat,indlon],axis=(1,2))
    for i in range(0,lstw):
        aux = np.nanmean(hs_gefs_hindcast[i,:,:,indlat,indlon],axis=1) # average ensemble members
        ind = np.where(aux>-999.)
        if np.size(ind)>int(np.floor(len(indt)/2)) :
            ind = ind[0]
            hs_gefs_hindcast_t[i] = np.nanmean(np.sort(aux[ind])[-nmax::])

        del ind,aux

    #  -- Probabilistic Forecast array --
    gspws=int(np.floor(wconfig['spws']/np.diff(latm).mean())/2)
    # U10
    prob_u10_gefs_forecast = np.zeros((lstw,len(qlev_wnd)),'f')*np.nan
    u10_gefs_forecast_t = np.zeros((lstw),'f')*np.nan
    for i in range(0,lstw):
        fmod=np.array(u10_gefs_forecast[i,:,:,(indlat-gspws):(indlat+gspws+1),:][:,:,:,(indlon-gspws):(indlon+gspws+1)])
        fmod=tsort(fmod)
        fmod=np.copy(fmod[-nmax::,:,:,:])
        fmod=np.array(fmod.reshape(nmax*len(ensm),fmod.shape[2],fmod.shape[3]))
        fmod=np.array(fmod.reshape(nmax*len(ensm),fmod.shape[1]*fmod.shape[2]))
        fmod=np.sort(fmod,axis=1)
        ind=np.where(np.mean(fmod,axis=0)>=0.)
        if np.size(ind)>0:
            fmod=np.array(fmod[:,ind[0]])            
            fmod=np.array(fmod[:,int(np.floor(fmod.shape[1]*(wconfig['spctl']/100)))::])
            fmod=fmod.reshape(fmod.shape[0]*fmod.shape[1])
            for j in range(0,len(qlev_wnd)):
                prob_u10_gefs_forecast[i,j] = np.size(fmod[fmod>qlev_wnd[j]]) / np.size(fmod[fmod>0.])
                if j==0:
                    u10_gefs_forecast_t[i] = np.nanmean(fmod)

        if i==0:
            fmod_result_u10 = np.array([fmod])
        else:
            fmod_result_u10 = np.append(fmod_result_u10,np.array([fmod]),axis=0)


    # Hs
    prob_hs_gefs_forecast = np.zeros((lstw,len(qlev_hs)),'f')*np.nan
    hs_gefs_forecast_t = np.zeros((lstw),'f')*np.nan
    for i in range(0,lstw):
        fmod=np.array(hs_gefs_forecast[i,:,:,(indlat-gspws):(indlat+gspws+1),:][:,:,:,(indlon-gspws):(indlon+gspws+1)])
        fmod = tsort(fmod)
        fmod=np.copy(fmod[-nmax::,:,:,:])
        fmod=np.array(fmod.reshape(nmax*len(ensm),fmod.shape[2],fmod.shape[3]))
        fmod=np.array(fmod.reshape(nmax*len(ensm),fmod.shape[1]*fmod.shape[2]))
        fmod=np.sort(fmod,axis=1)
        ind=np.where(np.mean(fmod,axis=0)>=0.)
        if np.size(ind)>0:
            fmod=np.array(fmod[:,ind[0]])            
            fmod=np.array(fmod[:,int(np.floor(fmod.shape[1]*(wconfig['spctl']/100)))::])
            fmod=fmod.reshape(fmod.shape[0]*fmod.shape[1])
            for j in range(0,len(qlev_hs)):
                prob_hs_gefs_forecast[i,j] = np.size(fmod[fmod>qlev_hs[j]]) / np.size(fmod[fmod>0.])
                if j==0:
                    hs_gefs_forecast_t[i] = np.nanmean(fmod)

        if i==0:
            fmod_result_hs = np.array([fmod])
        else:
            fmod_result_hs = np.append(fmod_result_hs,np.array([fmod]),axis=0)





    #  =====  =====  ===== 

    #  ++++++++ VALIDATION +++++++++++++++

    # ===== Categorical approach, verification of binary events =====

    # U10
    for i in range(0,len(qlev_wnd)):

        # Deterministic, matching positions and time
        if i==0:
            mop=ModelObsPlot(ru10_gdas[:,indlat,indlon],ru10_ndbc,ftag=ftag+"DeterministicVal_GDAS_U10_")
            mop.scatterplot(); mop.pdf(); mop.qqplot()
            mop=ModelObsPlot(np.nanmean(ru10_gefs_hindcast[:,:,indlat,indlon],axis=1),ru10_ndbc,ftag=ftag+"DeterministicVal_GEFShindcast_U10_")
            mop.scatterplot(); mop.pdf(); mop.qqplot()

        # Fuzzy verification
        # GDAS
        ceval_gdas = categorical_eval(u10_gdas_t,u10_ndbc_t,qlev_wnd[i])
        # GEFS hindcast
        ceval_gef_shindcast = categorical_eval(u10_gefs_hindcast_t,u10_ndbc_t,qlev_wnd[i])
        # GEFS forecast
        ceval_gef_forecast = categorical_eval(u10_gefs_forecast_t,u10_ndbc_t,qlev_wnd[i])

        dict_list = [
            {'Name': 'gdas', **ceval_gdas},
            {'Name': 'gefs_hindcast', **ceval_gef_shindcast},
            {'Name': 'gefs_forecast', **ceval_gef_forecast}
        ]

        # Convert the list of dictionaries into a pandas DataFrame
        df = pd.DataFrame(dict_list)
        df.to_csv(ftag+"CategoricalVal_InitialComparison_U10_Lev"+repr(np.round(fqlev_wnd[i],2))+".csv",sep='\t', index=False)

    # Hs
    for i in range(0,len(qlev_hs)):

        # Deterministic, matching positions and time
        if i==0:
            mop=ModelObsPlot(rhs_gdas[:,indlat,indlon],rhs_ndbc,ftag=ftag+"DeterministicVal_GDAS_Hs_")
            mop.scatterplot(); mop.pdf(); mop.qqplot()
            mop=ModelObsPlot(np.nanmean(rhs_gefs_hindcast[:,:,indlat,indlon],axis=1),rhs_ndbc,ftag=ftag+"DeterministicVal_GEFShindcast_Hs_")
            mop.scatterplot(); mop.pdf(); mop.qqplot()

        # Fuzzy verification
        # GDAS
        ceval_gdas = categorical_eval(hs_gdas_t,hs_ndbc_t,qlev_hs[i])
        # GEFS hindcast
        ceval_gef_shindcast = categorical_eval(hs_gefs_hindcast_t,hs_ndbc_t,qlev_hs[i])
        # GEFS forecast
        ceval_gef_forecast = categorical_eval(hs_gefs_forecast_t,hs_ndbc_t,qlev_hs[i])

        dict_list = [
            {'Name': 'gdas', **ceval_gdas},
            {'Name': 'gefs_hindcast', **ceval_gef_shindcast},
            {'Name': 'gefs_forecast', **ceval_gef_forecast}
        ]

        # Convert the list of dictionaries into a pandas DataFrame
        df = pd.DataFrame(dict_list)
        df.to_csv(ftag+"CategoricalVal_InitialComparison_Hs_Lev"+repr(int(fqlev_hs[i]))+".csv",sep='\t', index=False)


    # Binary Categorical array
    fprob_u10_gefs_forecast = np.zeros((lstw,len(qlev_wnd),lplev),'f')*np.nan
    fprob_hs_gefs_forecast = np.zeros((lstw,len(qlev_hs),lplev),'f')*np.nan
    for i in range(0,lstw):
        # U10
        for j in range(0,len(qlev_wnd)):
            for k in range(0,lplev):
                if prob_u10_gefs_forecast[i,j] >= plevels[k]:
                    fprob_u10_gefs_forecast[i,j,k] = int(1)
                else:
                    fprob_u10_gefs_forecast[i,j,k] = int(0)

        # Hs
        for j in range(0,len(qlev_hs)):
            for k in range(0,lplev):
                if prob_hs_gefs_forecast[i,j] >= plevels[k]:
                    fprob_hs_gefs_forecast[i,j,k] = int(1)
                else:
                    fprob_hs_gefs_forecast[i,j,k] = int(0)


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
    for i in range(0,len(qlev_wnd)):
        fftag=ftag+"ProbEvents_NDBC_U10_"+repr(np.round(fqlev_wnd[i],2))
        reliability_curve(prob_u10_gefs_forecast[:,i],u10_ndbc_t,qlev_wnd[i],fftag)

        fftag=ftag+"ProbEvents_GDAS_U10_"+repr(np.round(fqlev_wnd[i],2))
        reliability_curve(prob_u10_gefs_forecast[:,i],u10_gdas_t,qlev_wnd[i],fftag)

        fftag=ftag+"ProbEvents_GEFShindcast_U10_"+repr(np.round(fqlev_wnd[i],2))
        reliability_curve(prob_u10_gefs_forecast[:,i],u10_gefs_hindcast_t,qlev_wnd[i],fftag)

    # Hs
    for i in range(0,len(qlev_hs)):
        fftag=ftag+"ProbEvents_NDBC_Hs_"+repr(np.round(fqlev_hs[i],2))
        reliability_curve(prob_hs_gefs_forecast[:,i],hs_ndbc_t,qlev_hs[i],fftag)

        fftag=ftag+"ProbEvents_GDAS_Hs_"+repr(np.round(fqlev_hs[i],2))
        reliability_curve(prob_hs_gefs_forecast[:,i],hs_gdas_t,qlev_hs[i],fftag)

        fftag=ftag+"ProbEvents_GEFShindcast_Hs_"+repr(np.round(fqlev_hs[i],2))
        reliability_curve(prob_hs_gefs_forecast[:,i],hs_gefs_hindcast_t,qlev_hs[i],fftag)

    plt.close('all')

