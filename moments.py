# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 10:26:49 2019

This file comupte simulated moments + optionally 
plots some graphs

@author: Fabio
"""

import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.pyplot import plot, draw, show
import matplotlib.backends.backend_pdf
from numba import njit, vectorize

def moment(agents,draw=True):
#This function compute moments coming from the simulation
#Optionally it can also plot graphs about them. It is feeded with
#matrixes coming from simulations

    #Import simulated values
    assets_t=agents.setup.agrid_c[agents.iassets]
    iexo=agents.iexo
    state=agents.state
    theta_t=agents.setup.thetagrid_fine[agents.itheta]
    setup = agents.setup
    
   
    #As a first thing we unpack assets and theta
    N=len(state)
    
    #Get states codes
    state_codes = {name: i for i, name in enumerate(agents.setup.state_names)}
    
    ###########################################
    #Moments: FLS over time by Relationship
    ###########################################
    
    
    flsm=np.ones(agents.setup.pars['T'])
    flsc=np.ones(agents.setup.pars['T'])
    
    
    for t in range(agents.setup.pars['T']):
        
        pick = agents.state[:,t]==2        
        if pick.any(): flsm[t] = np.array(setup.ls_levels)[agents.ils_i[pick,t]].mean()
        pick = agents.state[:,t]==3
        if pick.any(): flsc[t] = np.array(setup.ls_levels)[agents.ils_i[pick,t]].mean()
        
        
    ###########################################
    #Moments: Variables over Age
    ###########################################
    
    relt=np.zeros((len(state_codes),agents.setup.pars['T']))
    ass_rel=np.zeros((len(state_codes),agents.setup.pars['T']))
    inc_rel=np.zeros((len(state_codes),agents.setup.pars['T']))
    
    
    
    print('relations')
    for ist,sname in enumerate(state_codes):
        for t in range(agents.setup.pars['T']):
            
            
            ftrend = agents.setup.pars['f_wage_trend'][t]
            mtrend = agents.setup.pars['m_wage_trend'][t]
            
            #Arrays for preparation
            is_state = (state[:,t]==ist)
            ind = np.where(is_state)[0]
            zf,zm,psi=agents.setup.all_indices(t,iexo[ind,t])[1:4]
            
            #Relationship over time
            relt[ist,t]=np.sum(is_state)
            
            #Assets over time   
            ass_rel[ist,t]=np.mean(assets_t[ind,t])
            
            #Income over time
            if sname=="Female, single":
                inc_rel[ist,t]=np.mean(np.exp(agents.setup.exogrid.zf_t[t][zf]  + ftrend ))
                
            elif sname=="Male, single":
                 inc_rel[ist,t]=np.mean(np.exp(agents.setup.exogrid.zf_t[t][zm] + mtrend))
                
            elif sname=="Couple, C" or sname=="Couple, M":
                 inc_rel[ist,t]=np.mean(np.exp(agents.setup.exogrid.zf_t[t][zf] + ftrend)+np.exp(agents.setup.exogrid.zf_t[t][zm] + mtrend))
    
            else:
            
               print('Error: No relationship chosen')
               
              
    print('spells new')


    nspells = (state[:,1:]!=state[:,:-1]).astype(np.int).sum(axis=1).max() + 1
    
    state_beg = -1*np.ones((N,nspells),dtype=np.int8)
    time_beg = -1*np.ones((N,nspells),dtype=np.bool)
    did_end = np.zeros((N,nspells),dtype=np.bool)
    state_end = -1*np.ones((N,nspells),dtype=np.int8)
    time_end = -1*np.ones((N,nspells),dtype=np.bool)
    sp_length = -1*np.ones((N,nspells),dtype=np.int16)
    is_spell = np.zeros((N,nspells),dtype=np.bool)
    
    state_beg[:,0] = 0 # THIS ASSUMES EVERYONE STARTS AS SINGLE
    time_beg[:,0] = 0 
    sp_length[:,0] = 1
    is_spell[:,0] = True
    ispell = np.zeros((N,),dtype=np.int8)
    
    for t in range(1,agents.setup.pars['T']):
        ichange = (state[:,t-1] != state[:,t])
        sp_length[~ichange,ispell[~ichange]] += 1
        
        if not np.any(ichange): continue
        
        did_end[ichange,ispell[ichange]] = True
        
        is_spell[ichange,ispell[ichange]+1] = True
        sp_length[ichange,ispell[ichange]+1] = 1 # if change then 1 year right
        state_end[ichange,ispell[ichange]] = state[ichange,t]
        time_end[ichange,ispell[ichange]] = t-1
        state_beg[ichange,ispell[ichange]+1] = state[ichange,t]  
        time_beg[ichange,ispell[ichange]] = t
        
        ispell[ichange] = ispell[ichange]+1
        
        
    allspells_beg = state_beg[is_spell]
    allspells_len = sp_length[is_spell]
    allspells_end = state_end[is_spell] # may be -1 if not ended
    
        
    ###########################################
    #Moments: Spells Creation
    ###########################################
    #spells_type=list()
    #spells_length=list()
    #spells_end=list()

    #@njit
    #def loop(spells_type,spells_length,spells_end):
    
#    print('spells old')
    
#    for n in range(N):
#        for ist in range(4):
#            for t in range(agents.setup.pars['T']-1):
#                
#                if t==0:
#                    leng=0 # I think we need to assign 1 in the first period
                            
#                 
#                if (leng>=0) and (state[n,t]==ist): 
#                    leng=leng+1
#                
#                if (leng>0 and state[n,t]!=ist) or (t==agents.setup.pars['T']-2): 
#                   # I did not understand this logic, please look again
#                   # as I understand we need to compare previous state with the 
#                   # current one, though I do not see how this is done    
#                   # (please elaborate, I might be wrong)
#    
#                    spells_type=[ist] + spells_type
#                    spells_length=[leng] + spells_length
#                    spells_end=[state[n,t]] + spells_end
#                    leng=0
#                    
#                    
#    
    
     #   return spells_type,spells_length,spells_end
            
    #spells_type,spells_length,spells_end=loop(spells_type,spells_length,spells_end)
    #spells_type.reverse()
    #spells_length.reverse()
    #spells_end.reverse()
    #spells=np.array([np.array(spells_type),np.array(spells_length),np.array(spells_end)]).T
    
    
    # TO FABIO: I did not understand the format of spells above. Please check.
    # I think in your format if the spell did not end the 3rd column is equal
    # to the first column, so I added this like
    allspells_end[allspells_end==-1] = allspells_beg[allspells_end==-1] # FIXME
    
    spells = np.stack((allspells_beg,allspells_len,allspells_end),axis=1)
    
    
    
    
    
    print('end spells')
   
    # TO FABIO: I did not edit things anything past this point.
    # Please look why my spells are different than yours. 
    
    
    #Now divide spells by relationship nature
    
    for ist,sname in enumerate(state_codes):
        s = sname.replace(',', '')
        s = s.replace(' ', '')
        
        
        is_state= (spells[:,0]==ist)
        
        if not is_state.any(): continue
            
        ind = np.where(is_state)[0]
        globals()['spells_t'+s]=spells[ind,:]
        is_state= (globals()['spells_t'+s][:,1]!=0)
        ind = np.where(is_state)[0]
        globals()['spells_'+s]=globals()['spells_t'+s][ind,:]
    
    
    ##################################
    # Construct the Hazard functions
    #################################
        
    #Hazard of Divorce
    hazd=list()
    lgh=len(spells_CoupleM[:,0])
    for t in range(agents.setup.pars['T']):
        
        cond=spells_CoupleM[:,1]==t+1
        temp=spells_CoupleM[cond,2]
        cond1=temp!=2
        temp1=temp[cond1]
        if lgh>0:
            haz1=len(temp1)/lgh
            lgh=lgh-len(temp)
        else:
            haz1=0.0
        hazd=[haz1]+hazd
        
    hazd.reverse()
    hazd=np.array(hazd).T 
    
    #Hazard of Separation
    hazs=list()
    lgh=len(spells_CoupleC[:,0])
    for t in range(agents.setup.pars['T']):
        
        cond=spells_CoupleC[:,1]==t+1
        temp=spells_CoupleC[cond,2]
        cond1=temp==0
        temp1=temp[cond1]
        if lgh>0:
            haz1=len(temp1)/lgh
            lgh=lgh-len(temp)
        else:
            haz1=0.0
        hazs=[haz1]+hazs
        
    hazs.reverse()
    hazs=np.array(hazs).T
    
    #Hazard of Marriage (Cohabitation spells)
    hazm=list()
    lgh=len(spells_CoupleC[:,0])
    for t in range(agents.setup.pars['T']):
        
        cond=spells_CoupleC[:,1]==t+1
        temp=spells_CoupleC[cond,2]
        cond1=temp==2
        temp1=temp[cond1]
        if lgh>0:
            haz1=len(temp1)/lgh
            lgh=lgh-len(temp)
        else:
            haz1=0.0
        hazm=[haz1]+hazm
        
    hazm.reverse()
    hazm=np.array(hazm).T
    
    #Singles: Marriage vs. cohabitation transition
    #spells_s=np.append(spells_Femalesingle,spells_Malesingle,axis=0)
    spells_s = spells_Femalesingle
    cond=spells_s[:,2]>1
    spells_sc=spells_s[cond,2]
    condm=spells_sc==2
    sharem=len(spells_sc[condm])/max(len(spells_sc),0.0001)
    
    if draw:
    
        #Print something useful for debug and rest
        print('The share of singles choosing marriage is {0:.2f}'.format(sharem))
        cond=(state<2)
        if assets_t[cond].size:
            print('The max level of assets for singles is {:.2f}, the grid upper bound is {:.2f}'.format(np.amax(assets_t[cond]),max(agents.setup.agrid_s)))
        cond=(state>1)
        if assets_t[cond].size:
            print('The max level of assets for couples is {:.2f}, the grid upper bound is {:.2f}'.format(np.amax(assets_t[cond]),max(agents.setup.agrid_c)))
        
        #Setup a file for the graphs
        pdf = matplotlib.backends.backend_pdf.PdfPages("moments_graphs.pdf")
        
        
        #############################################
        # Hazard of Divorce, Separation and Marriage
        #############################################
        fig = plt.figure()
        f1=fig.add_subplot(2,1,1)
        
        plt.plot(np.array(range(agents.setup.pars['T'])), hazd,markersize=6, label='Hazard of Divorce')
        plt.plot(np.array(range(agents.setup.pars['T'])), hazs,markersize=6, label='Hazard of Separation')
        plt.plot(np.array(range(agents.setup.pars['T'])), hazm,markersize=6, label='Hazard of Marriage')
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3),
                  fancybox=True, shadow=True, ncol=3, fontsize='x-small')
        #plt.legend(loc='upper left', shadow=True, fontsize='x-small')
        plt.xlabel('Duration')
        plt.ylabel('Hazard')
                   
        ##########################################
        # Assets Over the Live Cycle
        ########################################## 
        fig = plt.figure()
        f2=fig.add_subplot(2,1,1)
        
        for ist,sname in enumerate(state_codes):
            plt.plot(np.array(range(agents.setup.pars['T'])), ass_rel[ist,],color=print(ist/len(state_codes)),markersize=6, label=sname)
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3),
                  fancybox=True, shadow=True, ncol=len(state_codes), fontsize='x-small')
        #plt.legend(loc='upper left', shadow=True, fontsize='x-small')
        plt.xlabel('Time')
        plt.ylabel('Assets')
        
        ##########################################
        # Income Over the Live Cycle
        ########################################## 
        fig = plt.figure()
        f3=fig.add_subplot(2,1,1)
        
        for ist,sname in enumerate(state_codes):
          
            plt.plot(np.array(range(agents.setup.pars['T'])), inc_rel[ist,],color=print(ist/len(state_codes)),markersize=6, label=sname)
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3),
                  fancybox=True, shadow=True, ncol=len(state_codes), fontsize='x-small')
        plt.xlabel('Time')
        plt.ylabel('Income')
                
                
        ##########################################
        # Relationship Over the Live Cycle
        ##########################################       
        fig = plt.figure()
        f4=fig.add_subplot(2,1,1)
        for ist,sname in enumerate(state_codes):
            plt.plot([],[],color=print(ist/len(state_codes)), label=sname)
        plt.stackplot(np.array(range(agents.setup.pars['T'])),relt[0,]/N,relt[1,]/N,relt[2,]/N,relt[3,]/N,
                      colors = ['b','y','g','r'])            
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3),
                  fancybox=True, shadow=True, ncol=len(state_codes), fontsize='x-small')
        plt.xlabel('Time')
        plt.ylabel('Share')
        
        ##########################################
        # FLS Over the Live Cycle
        ##########################################       
        fig = plt.figure()
        f5=fig.add_subplot(2,1,1)

        plt.plot(np.array(range(agents.setup.pars['T'])), flsm,color='r', label='Marriage') 
        plt.plot(np.array(range(agents.setup.pars['T'])), flsc,color='k', label='Cohabitation')          
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3),
                  fancybox=True, shadow=True, ncol=len(state_codes), fontsize='x-small')
        plt.xlabel('Time')
        plt.ylabel('FLS')

        ##########################################
        # Put graphs together
        ########################################## 
        #show()
        for fig in range(1, plt.gcf().number + 1): ## will open an empty extra figure :(
            pdf.savefig( fig )
       
        pdf.close()
        matplotlib.pyplot.close("all")
        
        