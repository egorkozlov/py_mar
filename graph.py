# -*- coding: utf-8 -*-
"""
This file creates Graphs based on Policy and Value Functions

@author: Fabio
"""

import numpy as np
import dill as pickle
import matplotlib.pyplot as plt 
import gzip



def graphs(setup,ai,zfi,zmi,psii,ti,thi):
    # Import Value Funcrtion previously saved on File
    with gzip.open('name_model.pkl', 'rb') as file:
        Packed = pickle.load(file)
        
    ################################################
    # Unpack Stuff to Make it easier to Visualize
    ################################################
    T = setup.pars['T']
    agrid = setup.agrid
    agrids = setup.agrids
    zfg = setup.exogrid.zf_t[ti]
    zmg = setup.exogrid.zm_t[ti]
    psig = setup.exogrid.psi_t[ti]
   
    Vfs,cfs,sfs=np.empty([3,len(agrids), len(zfg),T])
    Vms,cms,sms=np.empty([3,len(agrids), len(zmg),T])
    Vfm,Vmm,cm,sm=np.empty([4,len(agrid), len(zfg),len(zmg),len(psig),T,setup.ntheta])
    Vfc,Vmc,cc,sc=np.empty([4,len(agrid), len(zfg),len(zmg),len(psig),T,setup.ntheta])
    
    #Single Women
    for t in range(T):
        for i in range(len(zfg)):
            for j in range(len(agrids)):
                Vfs[j,i,t]=Packed[t]['Female, single']['V'][j,i]
                cfs[j,i,t]=Packed[t]['Female, single']['c'][j,i]
                sfs[j,i,t]=Packed[t]['Female, single']['s'][j,i]
                
    #Single Men
    for t in range(T):
        for i in range(len(zmg)):
            for j in range(len(agrids)):
                Vms[j,i,t]=Packed[t]['Male, single']['V'][j,i]
                cms[j,i,t]=Packed[t]['Male, single']['c'][j,i]
                sms[j,i,t]=Packed[t]['Male, single']['s'][j,i]
                
    #Couples: Marriage+Cohabitation
    for t in range(T):
        for j in range(len(agrid)):
            for i in range(setup.pars['nexo']):
                
                #Get the indexes from zf,zm,psi
                zf,zm,psi=setup.all_indices(i)[1:4]
                
                #Marriage
                Vmm[j,zf,zm,psi,t,]=Packed[t]['Couple, M']['VM'][j,i,]
                Vfm[j,zf,zm,psi,t,]=Packed[t]['Couple, M']['VF'][j,i,]
                cm[j,zf,zm,psi,t,]=Packed[t]['Couple, M']['c'][j,i,]
                sm[j,zf,zm,psi,t,]=Packed[t]['Couple, M']['s'][j,i,]
                
                
                #Cohabitation
                Vmc[j,zf,zm,psi,t,]=Packed[t]['Couple, C']['VM'][j,i,]
                Vfc[j,zf,zm,psi,t,]=Packed[t]['Couple, C']['VF'][j,i,]
                cc[j,zf,zm,psi,t,]=Packed[t]['Couple, C']['c'][j,i,]
                sc[j,zf,zm,psi,t,]=Packed[t]['Couple, C']['s'][j,i,]
    
    
    #########################################
    # Additional Variables needed for graphs
    #########################################
    
    #Account for Single-Marriage and Single-Cohabit thresholds
    trem=np.array([100.0])
    trec=np.array([100.0])
    for i in range(len(psig)):
        if((Vfm[ai,zfi,zmi,i,ti,thi]-Vfs[ai,zfi,ti]>0.001) and (Vmm[ai,zfi,zmi,i,ti,thi]-Vms[ai,zmi,ti]>0.001) and trem>50.0):
            trem=psig[i]
        if((Vfc[ai,zfi,zmi,i,ti,thi]-Vfs[ai,zfi,ti]>0.001) and (Vmc[ai,zfi,zmi,i,ti,thi]-Vms[ai,zmi,ti]>0.001) and trec>50.0):
            trec=psig[i]
            
    tre=min(trem,trec)
    if tre>50.0:
        tre=max(psig)
    
    
    #####################################
    ################################
    ## Actually Construct graphs
    ################################
    ####################################
    
    
    ##########################################
    # Value Functions wrt Love
    ########################################## 
    fig = plt.figure()
    f1=fig.add_subplot(2,1,1)
    plt.plot(psig, Vmm[ai,zfi,zmi,0:len(psig),ti,thi],'bo',markersize=6, label='Man, Marriage')
    plt.plot(psig, Vmc[ai,zfi,zmi,0:len(psig),ti,thi],'b',linewidth=0.4, label='Man, Cohabitation')
    plt.plot(psig, Vfc[ai,zfi,zmi,0:len(psig),ti,thi],'r',linewidth=0.4, label='Women, Cohabitation')
    plt.plot(psig, Vfm[ai,zfi,zmi,0:len(psig),ti,thi],'r*',markersize=6,label='Women, Marriage')
    plt.axvline(x=tre, color='b', linestyle='--', label='Treshold Single-Couple')
    #plt.axvline(x=treb, color='b', linestyle='--', label='Tresh Bilateral')
    plt.xlabel('Love')
    plt.ylabel('Utility')
    #plt.title('Utility  Divorce costs: men=0.5, women=0.5')
    legend = plt.legend(loc='upper left', shadow=True, fontsize='x-small')
    
    ##########################################
    # Surplues of Marriage wrt Cohabitation, Love grid
    ##########################################
    
    #Generate Marriage Surplus wrt cohabitation + Value Functions
    surpM = [None] * len(psig)
    surpW = [None] * len(psig)
   
    for i in range(len(psig)):
        surpM[i]=max(Vmm[ai,zfi,zmi,i,ti,thi]-Vmc[ai,zfi,zmi,i,ti,thi],0.0)
        surpW[i]=max(Vfm[ai,zfi,zmi,i,ti,thi]-Vfc[ai,zfi,zmi,i,ti,thi],0.0)

    
    #Graph for the Surplus
    zero = np.array([0.0] * psig)
    fig2 = plt.figure()
    plt.plot(psig, zero,'k',linewidth=1)
    plt.plot(psig, surpM,'b',linewidth=1.5, label='Man')
    plt.plot(psig, surpW,'r',linewidth=1.5, label='Women')
    plt.axvline(x=tre, color='b', linestyle='--', label='Treshold Single-Couple')
    #plt.axvline(x=treb, color='b', label='Tresh Bilateral')
    #plt.axvline(x=treu, color='r', linestyle='--', label='Tresh Unilateral')
    plt.ylim(-0.1*max(max(surpM),max(surpW),max(surpM),max(surpW)),1.1*max(max(surpM),max(surpW)))
    legend = plt.legend(loc='upper left', shadow=True, fontsize='x-small')
    plt.xlabel('Love')
    plt.ylabel('Marriage Surplus wrt Cohab.')
    print(333,surpM,surpW)
    
    
    ##########################################
    # Value Function and Assets
    ########################################## 
    fig = plt.figure()
    f3=fig.add_subplot(2,1,1)
    plt.plot(agrid, Vmm[0:len(agrid),zfi,zmi,psii,ti,thi],'bo',markersize=6,markevery=5, label='Man, Marriage')
    plt.plot(agrid, Vmc[0:len(agrid),zfi,zmi,psii,ti,thi],'b',linewidth=0.4, label='Man, Cohabitation')
    plt.plot(agrid, Vfc[0:len(agrid),zfi,zmi,psii,ti,thi],'r',linewidth=0.4, label='Women, Cohabitation')
    plt.plot(agrid, Vfm[0:len(agrid),zfi,zmi,psii,ti,thi],'r*',markersize=6,markevery=5,label='Women, Marriage')
    #plt.axvline(x=treb, color='b', linestyle='--', label='Tresh Bilateral')
    plt.ylabel('Utility')
    plt.xlabel('Assets')
    #plt.title('Utility  Divorce costs: men=0.5, women=0.5')
    legend = plt.legend(loc='upper left', shadow=True, fontsize='x-small')
    
    ##########################################
    # Surplues of Marriage wrt Cohabitation, Asset
    ##########################################
    
    #Generate Marriage Surplus wrt cohabitation + Value Functions
    surpM = [None] * len(agrid)
    surpW = [None] * len(agrid)
   
    for i in range(len(agrid)):
        surpM[i]=max(Vmm[i,zfi,zmi,psii,ti,thi]-Vmc[i,zfi,zmi,psii,ti,thi],0.0)
        surpW[i]=max(Vfm[i,zfi,zmi,psii,ti,thi]-Vfc[i,zfi,zmi,psii,ti,thi],0.0)

    
    #Graph for the Surplus
    zero = np.array([0.0] * agrid)
    fig4 = plt.figure()
    plt.plot(agrid, zero,'k',linewidth=1)
    plt.plot(agrid, surpM,'b',linewidth=1.5, label='Man')
    plt.plot(agrid, surpW,'r',linewidth=1.5, label='Women')
    #plt.axvline(x=treb, color='b', label='Tresh Bilateral')
    #plt.axvline(x=treu, color='r', linestyle='--', label='Tresh Unilateral')
    plt.ylim(-0.1*max(max(surpM),max(surpW),max(surpM),max(surpW)),1.1*max(max(surpM),max(surpW)))
    legend = plt.legend(loc='upper left', shadow=True, fontsize='x-small')
    plt.xlabel('Assets')
    plt.ylabel('Marriage Surplus wrt Cohab.')
    
    ##########################################
    # Consumption and Assets
    ########################################## 
    fig = plt.figure()
    f5=fig.add_subplot(2,1,1)
    plt.plot(agrid, cm[0:len(agrid),zfi,zmi,psii,ti,thi],'ko',markersize=6,markevery=1, label='Marriage')
    plt.plot(agrid, cc[0:len(agrid),zfi,zmi,psii,ti,thi],'r*',markersize=6,markevery=1, label='Cohabitation')
    #plt.plot(agrid, sc[0:len(agrid),zfi,zmi,psii,ti,thi],'k',linewidth=2.0,linestyle='--', label='Cohabitation')
    plt.plot(agrids, cms[0:len(agrids),zmi,ti],'b',linewidth=2.0,label='Men, Single')
    plt.plot(agrids, cfs[0:len(agrids),zfi,ti],'r',linewidth=2.0,linestyle='--', label='Women, Single')
    #plt.axvline(x=treb, color='b', linestyle='--', label='Tresh Bilateral')
    plt.ylabel('Consumption')
    plt.xlabel('Assets')
    #plt.title('Utility  Divorce costs: men=0.5, women=0.5')
    legend = plt.legend(loc='upper left', shadow=True, fontsize='x-small')
    
    
    ##########################################
    # Savings and Assets
    ########################################## 
    fig = plt.figure()
    f6=fig.add_subplot(2,1,1)
    plt.plot(agrid, sm[0:len(agrid),zfi,zmi,psii,ti,thi],'ko',markersize=6,markevery=1, label='Marriage')
    plt.plot(agrid, sc[0:len(agrid),zfi,zmi,psii,ti,thi],'r*',markersize=6,markevery=1, label='Cohabitation')
    plt.plot(agrid, agrid,'k',linewidth=1.0,linestyle='--')
    plt.plot(agrids, sms[0:len(agrids),zmi,ti],'b',linewidth=2.0,label='Men, Single')
    plt.plot(agrids, sfs[0:len(agrids),zfi,ti],'r',linewidth=2.0,linestyle='--', label='Women, Single')
    #plt.axvline(x=treb, color='b', linestyle='--', label='Tresh Bilateral')
    plt.ylabel('Savings')
    plt.xlabel('Assets')
    #plt.title('Utility  Divorce costs: men=0.5, women=0.5')
    legend = plt.legend(loc='upper left', shadow=True, fontsize='x-small')
    #print(111,cms[0:len(agrid),zmi,ti])
    
    ##########################################
    # Value Function and Pareto Weights
    ########################################## 
    fig = plt.figure()
    f7=fig.add_subplot(2,1,1)
    plt.plot(setup.thetagrid, Vmm[ai,zfi,zmi,psii,ti,0:len(setup.thetagrid)],'bo',markersize=6, label='Man, Marriage')
    plt.plot(setup.thetagrid, Vmc[ai,zfi,zmi,psii,ti,0:len(setup.thetagrid)],'b',linewidth=0.4, label='Man, Cohabitation')
    plt.plot(setup.thetagrid, Vfc[ai,zfi,zmi,psii,ti,0:len(setup.thetagrid)],'r',linewidth=0.4, label='Women, Cohabitation')
    plt.plot(setup.thetagrid, Vfm[ai,zfi,zmi,psii,ti,0:len(setup.thetagrid)],'r*',markersize=6,label='Women, Marriage')
    #plt.axvline(x=treb, color='b', linestyle='--', label='Tresh Bilateral')
    plt.ylabel('Utility')
    plt.xlabel('Pareto Weight-Women')
    #plt.title('Utility  Divorce costs: men=0.5, women=0.5')
    legend = plt.legend(loc='upper left', shadow=True, fontsize='x-small')
    
    ###########################################
    # Value of Marriage-Cohabitation over time
    ###########################################
    #Generate Marriage Surplus wrt cohabitation + Value Functions
    surpM = [None] * T
    surpW = [None] * T
   
    for i in range(T):
        surpM[i]=max(Vmm[ai,zfi,zmi,psii,i,thi]-Vmc[ai,zfi,zmi,psii,i,thi],0.0)
        surpW[i]=max(Vfm[ai,zfi,zmi,psii,i,thi]-Vfc[ai,zfi,zmi,psii,i,thi],0.0)

    
    #Graph for the Surplus
    zero = np.array([0.0] * np.array(range(T)))
    fig8 = plt.figure()
    plt.plot(np.array(range(T)), zero,'k',linewidth=1)
    plt.plot(np.array(range(T)), surpM,'b',linewidth=1.5, label='Man')
    plt.plot(np.array(range(T)), surpW,'r',linewidth=1.5, label='Women')
    #plt.axvline(x=treb, color='b', label='Tresh Bilateral')
    #plt.axvline(x=treu, color='r', linestyle='--', label='Tresh Unilateral')
    plt.ylim(-0.1*max(max(surpM),max(surpW),max(surpM),max(surpW)),1.1*max(max(surpM),max(surpW)))
    legend = plt.legend(loc='upper left', shadow=True, fontsize='x-small')
    plt.xlabel('Time')
    plt.ylabel('Marriage Surplus wrt Cohab.')
    
    ##########################################
    # Consumption over the Life Cycle
    ########################################## 
    
       
    return Packed,cfs