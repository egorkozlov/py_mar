# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 10:26:49 2019
 
This file comupte simulated moments + optionally
plots some graphs
 
@author: Fabio
"""
 
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.pyplot import plot, draw, show
import matplotlib.backends.backend_pdf
import pickle
 
def moment(mdl,draw=True):
#This function compute moments coming from the simulation
#Optionally it can also plot graphs about them. It is feeded with
#matrixes coming from simulations
 
    agents = mdl.agents
    #Import simulated values
    assets_t_ifc=agents.setup.agrid_c[agents.iassets]
    assets_t_ifs=agents.setup.agrid_s[agents.iassets]
    
    
    iexo=agents.iexo
    state=agents.state
    theta_t=agents.setup.thetagrid_fine[agents.itheta]
    setup = agents.setup
    
    assets_t = np.zeros_like(assets_t_ifc)
    assets_t[(state==0) | (state==1)] = assets_t_ifs[(state==0) | (state==1)]
    assets_t[(state==2) | (state==3)] = assets_t_ifc[(state==2) | (state==3)]

        
    mdl.moments = dict()
     
     
    ##########################################
    #START COMPUTATION OF SIMULATED MOMENTS
    #########################################
    
    #As a first thing we unpack assets and theta
    N=len(state)
     
    #Get states codes
    state_codes = {name: i for i, name in enumerate(agents.setup.state_names)}
    
    
    ###########################################
    #Moments: Construction of Spells
    ###########################################
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
     
    # If the spell did not end mark it as ended with the state at its start
    allspells_end[allspells_end==-1] = allspells_beg[allspells_end==-1]
     
    spells = np.stack((allspells_beg,allspells_len,allspells_end),axis=1)
     
     
    #Now divide spells by relationship nature
    all_spells=dict()
    for ist,sname in enumerate(state_codes):
        #s = sname.replace(',', '')
        #s = s.replace(' ', '')
         
         
        is_state= (spells[:,0]==ist)
         
        #if not is_state.any(): continue
             
    
        all_spells[sname]=spells[is_state,:]

        is_state= (all_spells[sname][:,1]!=0)
        all_spells[sname]=all_spells[sname][is_state,:]
        
     
     
    ##################################
    # Construct the Hazard functions
    #################################
         
    #Hazard of Divorce
    hazd=list()
    lgh=len(all_spells['Couple, M'][:,0])
    for t in range(agents.setup.pars['Tret']):
         
        cond=all_spells['Couple, M'][:,1]==t+1
        temp=all_spells['Couple, M'][cond,2]
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
    lgh=len(all_spells['Couple, C'][:,0])
    for t in range(agents.setup.pars['Tret']):
         
        cond=all_spells['Couple, C'][:,1]==t+1
        temp=all_spells['Couple, C'][cond,2]
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
    lgh=len(all_spells['Couple, C'][:,0])
    for t in range(agents.setup.pars['Tret']):
         
        cond=all_spells['Couple, C'][:,1]==t+1
        temp=all_spells['Couple, C'][cond,2]
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
     
     
    mdl.moments['hazard sep'] = hazs
    mdl.moments['hazard div'] = hazd
    mdl.moments['hazard mar'] = hazm
     
     
 
     
    #Singles: Marriage vs. cohabitation transition
    #spells_s=np.append(spells_Femalesingle,spells_Malesingle,axis=0)
    spells_s =all_spells['Female, single']
    cond=spells_s[:,2]>1
    spells_sc=spells_s[cond,2]
    condm=spells_sc==2
    sharem=len(spells_sc[condm])/max(len(spells_sc),0.0001)
     
    ###########################################
    #Moments: FLS over time by Relationship
    ###########################################
     
     
    flsm=np.ones(agents.setup.pars['Tret'])
    flsc=np.ones(agents.setup.pars['Tret'])
     
     
    for t in range(agents.setup.pars['Tret']):
         
        pick = agents.state[:,t]==2       
        if pick.any(): flsm[t] = np.array(setup.ls_levels)[agents.ils_i[pick,t]].mean()
        pick = agents.state[:,t]==3
        if pick.any(): flsc[t] = np.array(setup.ls_levels)[agents.ils_i[pick,t]].mean()
         
     
         
    mdl.moments['flsm'] = flsm
    mdl.moments['flsc'] = flsc
    
    
    ###########################################
    #Sample selection
    ###########################################
    
    #Sample Selection to replicate the fact that
    #in NSFH wave two cohabitning couples were
    #excluded.
    #Birth cohorts: 45-55
    #Second wave of NLSFH:1992-1994.
    #
    #Assume that people are interviewd in 1993 and that age is uniformly
    #distributed. Clearly we can adjust this later on.
    
    
    
    #First cut the first two periods give new 'length'
    lenn=agents.setup.pars['T']-agents.setup.pars['Tbef']
    assets_t=assets_t[:,agents.setup.pars['Tbef']:agents.setup.pars['T']]
    iexo=iexo[:,agents.setup.pars['Tbef']:agents.setup.pars['T']]
    state=state[:,agents.setup.pars['Tbef']:agents.setup.pars['T']]
    theta_t=theta_t[:,agents.setup.pars['Tbef']:agents.setup.pars['T']]
    
    
    #Now drop observation to mimic the actual data gathering process
    keep=(assets_t[:,0]>-1)
    age_radius=int(10/agents.setup.pars['py'])

    for i in range(age_radius):
        keep[int(len(state[:,0])/age_radius*i):int(len(state[:,0])/age_radius*(i+1))]=\
        (state[int(len(state[:,0])/age_radius*i):int(len(state[:,0])/age_radius*(i+1)),int(28/agents.setup.pars['py']-i)]!=3)
    assets_t=assets_t[keep,]
    iexo=iexo[keep,]
    state=state[keep,]
    theta_t=theta_t[keep,]
     
    ###########################################
    #Moments: Variables over Age
    ###########################################
    
    #Update N to the new sample size
    N=len(state)
     
    relt=np.zeros((len(state_codes),lenn))
    relt1=np.zeros((len(state_codes),lenn))
    ass_rel=np.zeros((len(state_codes),lenn))
    inc_rel=np.zeros((len(state_codes),lenn))
     
     
     
    for ist,sname in enumerate(state_codes):
        for t in range(lenn):
             
            s=agents.setup.pars['Tbef']+t 
            ftrend = agents.setup.pars['f_wage_trend'][s]
            mtrend = agents.setup.pars['m_wage_trend'][s]
             
            #Arrays for preparation
            is_state = (np.any(state[:,0:t]==ist,1))       
            is_state1 = (state[:,t]==ist)
            if t<1:
                is_state=is_state1
            ind = np.where(is_state)[0]
            ind1 = np.where(is_state1)[0]
             
            if not (np.any(is_state) or np.any(is_state1)): continue
         
            zf,zm,psi=agents.setup.all_indices(t,iexo[ind1,t])[1:4]
             
            #Relationship over time
            relt[ist,t]=np.sum(is_state)
            relt1[ist,t]=np.sum(is_state1)
             
            #Assets over time  
            ass_rel[ist,t]=np.mean(assets_t[ind1,t])
            
             
            #Income over time
            if sname=="Female, single":
                inc_rel[ist,t]=np.mean(np.exp(agents.setup.exogrid.zf_t[s][zf]  + ftrend ))
                 
            elif sname=="Male, single":
                 inc_rel[ist,t]=np.mean(np.exp(agents.setup.exogrid.zf_t[s][zm] + mtrend))
                 
            elif sname=="Couple, C" or sname=="Couple, M":
                 inc_rel[ist,t]=np.mean(np.exp(agents.setup.exogrid.zf_t[s][zf] + ftrend)+np.exp(agents.setup.exogrid.zf_t[s][zm] + mtrend))
     
            else:
             
               print('Error: No relationship chosen')
              
    #Now, before saving the moments, take interval of 5 years
    # if (agents.setup.pars['Tret']>=agents.setup.pars['Tret']):        
    reltt=relt[:,0:agents.setup.pars['Tret']-agents.setup.pars['Tbef']+1]
    years=np.linspace(20,50,7)
    years_model=np.linspace(20,50,30/agents.setup.pars['py'])
    
    #Find the right entries for creating moments
    pos=list()
    for j in range(len(years)):
        pos=pos+[np.argmin(np.abs(years_model-years[j]))]
    
    #Approximation if more than 5 years in one period
    if len(pos)<7:
        for i in range(7-len(pos)):
            pos=pos+[pos[-1]]
    pos=np.array(pos)
    
    
    
    reltt=reltt[:,pos]
    #else:
     #   reltt=relt
        
    mdl.moments['share single'] = reltt[0,:]/N
    mdl.moments['share mar'] = reltt[2,:]/N
    mdl.moments['share coh'] = reltt[3,:]/N
               
     
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
         
        #################
        #Get data moments
        #################
         
        #Get Data Moments
        with open('moments.pkl', 'rb') as file:
            packed_data=pickle.load(file)
         
            #Unpack Moments (see data_moments.py to check if changes)
            #(hazm,hazs,hazd,mar,coh,fls_ratio,W)
            hazm_d=packed_data[0]
            hazs_d=packed_data[1]
            hazd_d=packed_data[2]
            mar_d=packed_data[3]
            coh_d=packed_data[4]
            fls_d=np.ones(1)*packed_data[5]
            hazm_i=packed_data[7]
            hazs_i=packed_data[8]
            hazd_i=packed_data[9]
            mar_i=packed_data[10]
            coh_i=packed_data[11]
            fls_i=np.ones(1)*packed_data[12]
 
         
         
        #############################################
        # Hazard of Divorce
        #############################################
        fig = plt.figure()
        f1=fig.add_subplot(2,1,1)
        lg=min(len(hazd_d),len(hazd))
        if lg<2:
            one='o'
            two='o'
        else:
            one='r'
            two='b'
        plt.plot(np.array(range(lg)), hazd[0:lg],one, linestyle='--',linewidth=1.5, label='Hazard of Divorce - S')
        plt.plot(np.array(range(lg)), hazd_d[0:lg],two,linewidth=1.5, label='Hazard of Divorce - D')
        plt.fill_between(np.array(range(lg)), hazd_i[0,0:lg], hazd_i[1,0:lg],alpha=0.2,facecolor='b')
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3),
                  fancybox=True, shadow=True, ncol=3, fontsize='x-small')
        plt.ylim(ymin=0)
        #plt.legend(loc='upper left', shadow=True, fontsize='x-small')
        plt.xlabel('Duration')
        plt.ylabel('Hazard')
         
        #############################################
        # Hazard of Separation
        #############################################
        fig = plt.figure()
        f1=fig.add_subplot(2,1,1)
        lg=min(len(hazs_d),len(hazs))
        plt.plot(np.array(range(lg)), hazs[0:lg],one, linestyle='--',linewidth=1.5, label='Hazard of Separation - S')
        plt.plot(np.array(range(lg)), hazs_d[0:lg],two,linewidth=1.5, label='Hazard of Separation - D')
        plt.fill_between(np.array(range(lg)), hazs_i[0,0:lg], hazs_i[1,0:lg],alpha=0.2,facecolor='b')
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3),
                  fancybox=True, shadow=True, ncol=3, fontsize='x-small')
        plt.ylim(ymin=0)
        #plt.legend(loc='upper left', shadow=True, fontsize='x-small')
        plt.xlabel('Duration')
        plt.ylabel('Hazard')
        
         
        #############################################
        # Hazard of Marriage
        #############################################
        fig = plt.figure()
        f1=fig.add_subplot(2,1,1)
        lg=min(len(hazm_d),len(hazm))

        plt.plot(np.array(range(lg)), hazm[0:lg],one, linestyle='--',linewidth=1.5, label='Hazard of Marriage - S')
        plt.plot(np.array(range(lg)), hazm_d[0:lg],two,linewidth=1.5, label='Hazard of Marriage - D')
        plt.fill_between(np.array(range(lg)), hazm_i[0,0:lg], hazm_i[1,0:lg],alpha=0.2,facecolor='b')
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3),
                  fancybox=True, shadow=True, ncol=3, fontsize='x-small')
        plt.ylim(ymin=0)
        #plt.legend(loc='upper left', shadow=True, fontsize='x-small')
        plt.xlabel('Duration')
        plt.ylabel('Hazard')
        
        ##########################################
        # Assets Over the Live Cycle
        ##########################################
        fig = plt.figure()
        f2=fig.add_subplot(2,1,1)
         
        for ist,sname in enumerate(state_codes):
            plt.plot(np.array(range(lenn)), ass_rel[ist,],color=print(ist/len(state_codes)),markersize=6, label=sname)
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
           
            plt.plot(np.array(range(lenn)), inc_rel[ist,],color=print(ist/len(state_codes)),markersize=6, label=sname)
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
        plt.stackplot(np.array(range(len(relt1[0,]))),relt1[0,]/N,relt1[1,]/N,relt1[2,]/N,relt1[3,]/N,
                      colors = ['b','y','g','r'])           
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3),
                  fancybox=True, shadow=True, ncol=len(state_codes), fontsize='x-small')
        plt.xlabel('Time')
        plt.ylabel('Share')
         
        ##########################################
        # Relationship and Data
        ##########################################      
        fig = plt.figure()
        f4=fig.add_subplot(2,1,1)
        lg=min(len(mar_d),len(relt[1,:]))
        plt.plot(np.array(range(lg)), mar_d[0:lg],'g',linewidth=1.5, label='Share Married - D')
        plt.fill_between(np.array(range(lg)), mar_i[0,0:lg], mar_i[1,0:lg],alpha=0.2,facecolor='g')
        plt.plot(np.array(range(lg)), reltt[2,0:lg]/N,'g',linestyle='--',linewidth=1.5, label='Share Married - S')
        plt.plot(np.array(range(lg)), coh_d[0:lg],'r',linewidth=1.5, label='Share Cohabiting - D')
        plt.fill_between(np.array(range(lg)), coh_i[0,0:lg], coh_i[1,0:lg],alpha=0.2,facecolor='r')
        plt.plot(np.array(range(lg)), reltt[3,0:lg]/N,'r',linestyle='--',linewidth=1.5, label='Share Cohabiting - S')
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3),
                  fancybox=True, shadow=True, ncol=len(state_codes), fontsize='x-small')
        plt.ylim(ymax=1.0)
        plt.xlabel('Time')
        plt.ylabel('Share')
         
        ##########################################
        # FLS Over the Live Cycle
        ##########################################      
        fig = plt.figure()
        f5=fig.add_subplot(2,1,1)
 
        plt.plot(np.array(range(agents.setup.pars['Tret'])), flsm,color='r', label='Marriage')
        plt.plot(np.array(range(agents.setup.pars['Tret'])), flsc,color='k', label='Cohabitation')         
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
         
        
