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
import pandas as pd  
import statsmodels.api as sm  
import statsmodels.formula.api as smf  
 
   
def moment(mdl,agents,draw=True,validation=False):  
#This function compute moments coming from the simulation  
#Optionally it can also plot graphs about them. It is feeded with  
#matrixes coming from simulations  
   
  
    #Import simulated values  
    assets_t=mdl.setup.agrid_c[agents.iassets] # FIXME  
    iexo=agents.iexo  
    state=agents.state  
    theta_t=mdl.setup.thetagrid_fine[agents.itheta]  
    setup = mdl.setup  
      
    
          
          
    moments = dict()  
       
       
    ##########################################  
    #START COMPUTATION OF SIMULATED MOMENTS  
    #########################################  
     
       
    #Create a file with the age of the change foreach person  
    changep=agents.policy_ind 
      
       
    #Get states codes  
    state_codes = {name: i for i, name in enumerate(mdl.setup.state_names)}  
     
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
    assets_t=assets_t[:,0:mdl.setup.pars['T']]  
    iexo=iexo[:,0:mdl.setup.pars['T']]  
    state=state[:,0:mdl.setup.pars['T']]  
    theta_t=theta_t[:,0:mdl.setup.pars['T']]  
     
      
    ####################################################################  
    #Now drop observation to mimic the actual data gathering process  
    ####################################################################  
      
    #Get distribution of age conditional on cohabiting on the second wave  
    with open('age_sw.pkl', 'rb') as file:  
        age_sw=pickle.load(file)  
          
    keep=(assets_t[:,0]>-1)  
     
  
    summa=0.0  
    summa1=0.0  
    for i in age_sw:  
        summa+=age_sw[i]  
        keep[int(summa1*len(state[:,0])/sum(age_sw.values())):int(summa*len(state[:,0])/sum(age_sw.values()))]=(state[int(summa1*len(state[:,0])/sum(age_sw.values())):int(summa*len(state[:,0])/sum(age_sw.values())),int((i-20)/mdl.setup.pars['py'])]!=3)  
        
        summa1+=age_sw[i]  
    assets_t=assets_t[keep,]  
    iexo=iexo[keep,]  
    state=state[keep,]  
    theta_t=theta_t[keep,]  
    changep=changep[keep,] 
     
    index=np.array(np.linspace(1,len(state[:,0]),len(state[:,0]))-1,dtype=np.int16) 
     
    N=len(iexo[:,0]) 
      
    ################################################################### 
    #Get age we stop observing spells 
    ################################################################### 
    with open('age_sint.pkl', 'rb') as file:  
        age_sint=pickle.load(file)  
          
    aged=np.ones((state.shape)) 
     
  
    summa=0.0 
    summa1=0.0 
    for i in age_sint: 
        summa+=age_sint[int(i)] 
        aged[int(summa1*len(aged[:])/sum(age_sint.values())):int(summa*len(aged[:])/sum(age_sint.values()))]=round((i-20)/mdl.setup.pars['py'],0) 
        summa1+=age_sint[int(i)] 
     
    aged=np.array(aged,dtype=np.int16) 
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
    sp_person = -1*np.ones((N,nspells),dtype=np.int16)  
    is_unid = -1*np.ones((N,nspells),dtype=np.int16)  
    is_unid_end = -1*np.ones((N,nspells),dtype=np.int16)  
    is_unid_lim = -1*np.ones((N,nspells),dtype=np.int16)  
    n_spell = -1*np.ones((N,nspells),dtype=np.int16)  
    is_spell = np.zeros((N,nspells),dtype=np.bool)  
     
       
    state_beg[:,0] = 0 # THIS ASSUMES EVERYONE STARTS AS SINGLE  
    time_beg[:,0] = 0  
    sp_length[:,0] = 1  
    is_spell[:,0] = True  
    ispell = np.zeros((N,),dtype=np.int8)  
       
    for t in range(1,mdl.setup.pars['T']):  
        ichange = ((state[:,t-1] != state[:,t]))  
        sp_length[((~ichange)),ispell[((~ichange))]] += 1  
        #ichange = ((state[:,t-1] != state[:,t]) & (t<=aged[:,t]))  
        #sp_length[((~ichange) & (t<=aged[:,t])),ispell[((~ichange) & (t<=aged[:,t]))]] += 1  
           
        if not np.any(ichange): continue  
           
        did_end[ichange,ispell[ichange]] = True  
           
        is_spell[ichange,ispell[ichange]+1] = True  
        sp_length[ichange,ispell[ichange]+1] = 1 # if change then 1 year right  
        state_end[ichange,ispell[ichange]] = state[ichange,t]  
        sp_person[ichange,ispell[ichange]] = index[ichange] 
        time_end[ichange,ispell[ichange]] = t-1  
        state_beg[ichange,ispell[ichange]+1] = state[ichange,t]   
        time_beg[ichange,ispell[ichange]+1] = t  
        n_spell[ichange,ispell[ichange]+1]=ispell[ichange]+1 
        is_unid[ichange,ispell[ichange]+1]=changep[ichange,t] 
        is_unid_lim[ichange,ispell[ichange]+1]=changep[ichange,aged[ichange,0]] 
        is_unid_end[ichange,ispell[ichange]]=changep[ichange,t-1] 
         
           
        ispell[ichange] = ispell[ichange]+1  
           
           
    allspells_beg = state_beg[is_spell]  
    allspells_len = sp_length[is_spell]  
    allspells_end = state_end[is_spell] # may be -1 if not ended  
    allspells_timeb = time_beg[is_spell] 
    allspells_isunid=is_unid[is_spell] 
    allspells_isunidend=is_unid_end[is_spell] 
    allspells_isunidlim=is_unid_lim[is_spell] 
    allspells_person=sp_person[is_spell] 
    allspells_nspells=n_spell[is_spell] 
     
     
       
    # If the spell did not end mark it as ended with the state at its start  
    allspells_end[allspells_end==-1] = allspells_beg[allspells_end==-1]  
    allspells_isunidend[allspells_isunidend==-1] = allspells_isunidlim[allspells_isunidend==-1] 
    allspells_nspells[allspells_nspells==-1]=0 
    allspells_nspells=allspells_nspells+1 
      
    #Use this to construct hazards 
    spells = np.stack((allspells_beg,allspells_len,allspells_end),axis=1)  
     
    #Use this for empirical analysis 
    spells_empirical=np.stack((allspells_beg,allspells_timeb,allspells_len,allspells_end,allspells_nspells,allspells_isunid,allspells_isunidend),axis=1) 
    is_coh=((spells_empirical[:,0]==3) & (spells_empirical[:,5]==spells_empirical[:,6])) 
    spells_empirical=spells_empirical[is_coh,1:6] 
     
    
       
       
    #Now divide spells by relationship nature  
    all_spells=dict()  
    for ist,sname in enumerate(state_codes):  
  
        is_state= (spells[:,0]==ist)  
           
   
        all_spells[sname]=spells[is_state,:]  
  
        is_state= (all_spells[sname][:,1]!=0)  
        all_spells[sname]=all_spells[sname][is_state,:]  
          
          
    ############################################  
    #Construct sample of first relationships  
    ############################################  
 
      
    #Now define variables  
    rel_end = -1*np.ones((N,99),dtype=np.int16)  
    rel_age= -1*np.ones((N,99),dtype=np.int16)  
    rel_unid= -1*np.ones((N,99),dtype=np.int16)  
    rel_number= -1*np.ones((N,99),dtype=np.int16)  
    isrel = np.zeros((N,),dtype=np.int8)  
      
    for t in range(1,mdl.setup.pars['Tret']):  
          
        irchange = ((state[:,t-1] != state[:,t]) & ((state[:,t-1]==0) | (state[:,t-1]==1)))  
          
        if not np.any(ichange): continue  
      
        rel_end[irchange,isrel[irchange]]=state[irchange,t]  
        rel_age[irchange,isrel[irchange]]=t  
        rel_unid[irchange,isrel[irchange]]=changep[irchange,t]  
        rel_number[irchange,isrel[irchange]]=isrel[irchange]+1  
          
        isrel[irchange] = isrel[irchange]+1  
      
    #Get the final Variables  
    allrel_end=rel_end[(rel_end!=-1)]  
    allrel_age=rel_age[(rel_age!=-1)]  
    allrel_uni=rel_unid[(rel_unid!=-1)]  
    allrel_number=rel_number[(rel_number!=-1)]  
      
    #Get whetehr marraige  
    allrel_mar=np.zeros((allrel_end.shape))  
    allrel_mar[(allrel_end==2)]=1  
      
    #Create a Pandas Dataframe  
    data_rel=np.array(np.stack((allrel_mar,allrel_age,allrel_uni,allrel_number),axis=0).T,dtype=np.float64)  
    data_rel_panda=pd.DataFrame(data=data_rel,columns=['mar','age','uni','rnumber'])  
                     
  
      
    #Regression  
    try:  
        FE_ols = smf.ols(formula='mar ~ uni+C(rnumber)+C(age)', data = data_rel_panda.dropna()).fit()  
        beta_unid_s=FE_ols.params['uni']  
    except:  
        print('No data for unilateral divorce regression...')  
        beta_unid_s=0.0  
      
      
    moments['beta unid']=beta_unid_s   
     
    ################################################### 
    # Second regression for the length of cohabitation 
    ################################################### 
    data_coh_panda=pd.DataFrame(data=spells_empirical,columns=['age','duration','end','rel','uni'])  
     
    #Regression  
    try:  
    #FE_ols = smf.ols(formula='duration ~ uni+C(age)', data = data_coh_panda.dropna()).fit()  
    #beta_dur_s=FE_ols.params['uni']  
     
        from lifelines import CoxPHFitter 
        cph = CoxPHFitter() 
        data_coh_panda['age2']=data_coh_panda['age']**2 
        data_coh_panda['age3']=data_coh_panda['age']**3 
        data_coh_panda['rel2']=data_coh_panda['rel']**2 
        data_coh_panda['rel3']=data_coh_panda['rel']**3 
        #data_coh_panda=pd.get_dummies(data_coh_panda, columns=['age']) 
         
        #Standard Cox 
        data_coh_panda['endd']=1.0 
        data_coh_panda.loc[data_coh_panda['end']==3.0,'endd']=0.0 
        data_coh_panda1=data_coh_panda.drop(['end'], axis=1) 
        cox_join=cph.fit(data_coh_panda1, duration_col='duration', event_col='endd') 
        haz_join=cox_join.hazard_ratios_['uni'] 
         
        #Cox where risk is marriage 
        data_coh_panda['endd']=0.0 
        data_coh_panda.loc[data_coh_panda['end']==2.0,'endd']=1.0 
        data_coh_panda2=data_coh_panda.drop(['end'], axis=1) 
        cox_mar=cph.fit(data_coh_panda2, duration_col='duration', event_col='endd') 
        haz_mar=cox_mar.hazard_ratios_['uni'] 
         
        #Cox where risk is separatio 
        data_coh_panda['endd']=0.0 
        data_coh_panda.loc[data_coh_panda['end']==0.0,'endd']=1.0 
        data_coh_panda3=data_coh_panda.drop(['end'], axis=1) 
        cox_sep=cph.fit(data_coh_panda3, duration_col='duration', event_col='endd') 
        haz_sep=cox_sep.hazard_ratios_['uni'] 
         
    except:  
        print('No data for unilateral divorce regression...')  
        haz_sep=1.0
        haz_join=1.0
        haz_mar=1.0
       
    ##################################  
    # Construct the Hazard functions  
    #################################  
           
    #Hazard of Divorce  
    hazd=list()  
    lgh=len(all_spells['Couple, M'][:,0])  
    for t in range(mdl.setup.pars['T']):  
           
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
    for t in range(mdl.setup.pars['T']):  
           
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
    for t in range(mdl.setup.pars['T']):  
           
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
       
       
    moments['hazard sep'] = hazs  
    moments['hazard div'] = hazd  
    moments['hazard mar'] = hazm  
      
  
   
       
    #Singles: Marriage vs. cohabitation transition  
    #spells_s=np.append(spells_Femalesingle,spells_Malesingle,axis=0)  
    spells_s =all_spells['Female, single']  
    cond=spells_s[:,2]>1  
    spells_sc=spells_s[cond,2]  
    condm=spells_sc==2  
    sharem=len(spells_sc[condm])/max(len(spells_sc),0.0001)  
     
     
    #Cut the first two periods give new 'length'  
    lenn=mdl.setup.pars['T']-mdl.setup.pars['Tbef']  
    assets_t=assets_t[:,mdl.setup.pars['Tbef']:mdl.setup.pars['T']]  
    iexo=iexo[:,mdl.setup.pars['Tbef']:mdl.setup.pars['T']]  
    state=state[:,mdl.setup.pars['Tbef']:mdl.setup.pars['T']]  
    theta_t=theta_t[:,mdl.setup.pars['Tbef']:mdl.setup.pars['T']]  
       
    ###########################################  
    #Moments: FLS over time by Relationship  
    ###########################################  
       
       
    flsm=np.ones(mdl.setup.pars['Tret'])  
    flsc=np.ones(mdl.setup.pars['Tret'])  
       
       
    for t in range(mdl.setup.pars['Tret']):  
           
        pick = agents.state[:,t]==2         
        if pick.any(): flsm[t] = np.array(setup.ls_levels)[agents.ils_i[pick,t]].mean()  
        pick = agents.state[:,t]==3  
        if pick.any(): flsc[t] = np.array(setup.ls_levels)[agents.ils_i[pick,t]].mean()  
           
       
           
    moments['flsm'] = flsm  
    moments['flsc'] = flsc  
      
      
    
       
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
               
            s=mdl.setup.pars['Tbef']+t   
            ftrend = mdl.setup.pars['f_wage_trend'][s]  
            mtrend = mdl.setup.pars['m_wage_trend'][s]  
               
            #Arrays for preparation  
            is_state = (np.any(state[:,0:t]==ist,1))         
            is_state1 = (state[:,t]==ist)  
            if t<1:  
                is_state=is_state1  
            ind = np.where(is_state)[0]  
            ind1 = np.where(is_state1)[0]  
               
            if not (np.any(is_state) or np.any(is_state1)): continue  
           
            zf,zm,psi=mdl.setup.all_indices(t,iexo[ind1,t])[1:4]  
               
            #Relationship over time  
            relt[ist,t]=np.sum(is_state)  
            relt1[ist,t]=np.sum(is_state1)  
               
            #Assets over time    
            ass_rel[ist,t]=np.mean(assets_t[ind1,t])  
              
               
            #Income over time  
            if sname=="Female, single":  
                inc_rel[ist,t]=np.mean(np.exp(mdl.setup.exogrid.zf_t[s][zf]  + ftrend ))  
                   
            elif sname=="Male, single":  
                 inc_rel[ist,t]=np.mean(np.exp(mdl.setup.exogrid.zf_t[s][zm] + mtrend))  
                   
            elif sname=="Couple, C" or sname=="Couple, M":  
                 inc_rel[ist,t]=np.mean(np.exp(mdl.setup.exogrid.zf_t[s][zf] + ftrend)+np.exp(mdl.setup.exogrid.zf_t[s][zm] + mtrend))  
       
            else:  
               
               print('Error: No relationship chosen')  
                
    #Now, before saving the moments, take interval of 5 years  
    # if (mdl.setup.pars['Tret']>=mdl.setup.pars['Tret']):          
    reltt=relt[:,0:mdl.setup.pars['Tret']-mdl.setup.pars['Tbef']+1]  
    years=np.linspace(20,50,7)  
    years_model=np.linspace(20,50,30/mdl.setup.pars['py'])  
      
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
          
    moments['share single'] = reltt[0,:]/N  
    moments['share mar'] = reltt[2,:]/N  
    moments['share coh'] = reltt[3,:]/N  
                 
       
    if draw:  
       
        #Print something useful for debug and rest  
        print('The share of singles choosing marriage is {0:.2f}'.format(sharem))  
        cond=(state<2)  
        if assets_t[cond].size:  
            print('The max level of assets for singles is {:.2f}, the grid upper bound is {:.2f}'.format(np.amax(assets_t[cond]),max(mdl.setup.agrid_s)))  
        cond=(state>1)  
        if assets_t[cond].size:  
            print('The max level of assets for couples is {:.2f}, the grid upper bound is {:.2f}'.format(np.amax(assets_t[cond]),max(mdl.setup.agrid_c)))  
           
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
            hazm_d=packed_data['hazm']  
            hazs_d=packed_data['hazs']  
            hazd_d=packed_data['hazd']  
            mar_d=packed_data['emar']  
            coh_d=packed_data['ecoh']  
            fls_d=np.ones(1)*packed_data['fls_ratio']  
            beta_unid_d=np.ones(1)*packed_data['beta_unid']  
            hazm_i=packed_data['hazmi']  
            hazs_i=packed_data['hazsi']  
            hazd_i=packed_data['hazdi']  
            mar_i=packed_data['emari']  
            coh_i=packed_data['ecohi']  
            fls_i=np.ones(1)*packed_data['fls_ratioi']  
            beta_unid_i=np.ones(1)*packed_data['beta_unidi']  
   
           
           
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
        xa=(mdl.setup.pars['py']*np.array(range(len(relt1[0,])))+20) 
        for ist,sname in enumerate(state_codes):  
            plt.plot([],[],color=print(ist/len(state_codes)), label=sname)  
        plt.stackplot(xa,relt1[0,]/N,relt1[1,]/N,relt1[2,]/N,relt1[3,]/N,  
                      colors = ['b','y','g','r'])             
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3),  
                  fancybox=True, shadow=True, ncol=len(state_codes), fontsize='x-small')  
        plt.xlabel('Age')  
        plt.ylabel('Share')  
           
        ##########################################  
        # Relationship and Data  
        ##########################################        
        fig = plt.figure()  
        f4=fig.add_subplot(2,1,1)  
        lg=min(len(mar_d),len(relt[1,:]))  
        xa=(5*np.array(range(lg))+20) 
        plt.plot(xa, mar_d[0:lg],'g',linewidth=1.5, label='Share Married - D')  
        plt.fill_between(xa, mar_i[0,0:lg], mar_i[1,0:lg],alpha=0.2,facecolor='g')  
        plt.plot(xa, reltt[2,0:lg]/N,'g',linestyle='--',linewidth=1.5, label='Share Married - S')  
        plt.plot(xa, coh_d[0:lg],'r',linewidth=1.5, label='Share Cohabiting - D')  
        plt.fill_between(xa, coh_i[0,0:lg], coh_i[1,0:lg],alpha=0.2,facecolor='r')  
        plt.plot(xa, reltt[3,0:lg]/N,'r',linestyle='--',linewidth=1.5, label='Share Cohabiting - S')  
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3),  
                  fancybox=True, shadow=True, ncol=len(state_codes), fontsize='x-small')  
        plt.ylim(ymax=1.0)  
        plt.xlabel('Age')  
        plt.ylabel('Share')  
           
        ##########################################  
        # FLS Over the Live Cycle  
        ##########################################        
        fig = plt.figure()  
        f5=fig.add_subplot(2,1,1)  
        xa=(mdl.setup.pars['py']*np.array(range(mdl.setup.pars['Tret']))+20) 
        plt.plot(xa, flsm,color='r', label='Marriage')  
        plt.plot(xa, flsc,color='k', label='Cohabitation')           
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3),  
                  fancybox=True, shadow=True, ncol=len(state_codes), fontsize='x-small')  
        plt.xlabel('Age')  
        plt.ylabel('FLS')  
          
        ##########################################  
        # Histogram of Effect of Unilateral Divorce  
        ##########################################   
        fig = plt.figure()  
        f6=fig.add_subplot(2,1,1)  
           
          
        # create plot  
        x=np.array([0.25,0.75]) 
        y=np.array([beta_unid_d,beta_unid_s])  
        yerr=np.array([(beta_unid_i[1]-beta_unid_i[0])/2.0,0.0])  
        plt.axhline(linewidth=0.1, color='r')  
        plt.errorbar(x, y, yerr=yerr, fmt='o', elinewidth=0.03)  
        plt.ylabel('OLS Coefficient - UniD')  
        plt.xticks(x, ["Data","Simulation"] ) 
        plt.ylim(ymax=0.1)  
        plt.xlim(xmax=1.0,xmin=0.0)  
        #plt.xticks(index , ('Unilateral', 'Bilateral'))  
          
         
        ########################################################## 
        # Histogram of Unilateral Divorce on Cohabitation Length 
        ############################################################# 
        fig = plt.figure()  
        f6=fig.add_subplot(2,1,1)  
           
          
        # create plot  
        x=np.array([0.2,0.5,0.8]) 
        y=np.array([haz_join,haz_mar,haz_sep])  
        yerr=y*0.0  
        plt.axhline(y=1.0,linewidth=0.1, color='r')  
        plt.errorbar(x, y, yerr=yerr, fmt='o', elinewidth=0.03)  
        plt.ylabel('Relative Hazard - UniD vs. Bil')  
        plt.xticks(x, ["Overall Risk","Risk of Marriage","Risk of Separation"] ) 
        #plt.ylim(ymax=1.2,ymin=0.7)  
        plt.xlim(xmax=1.0,xmin=0.0)  
        #plt.xticks(index , ('Unilateral', 'Bilateral'))  
   
        ##########################################  
        # Put graphs together  
        ##########################################  
        #show()  
        for fig in range(1, plt.gcf().number + 1): ## will open an empty extra figure :(  
            pdf.savefig( fig )  
          
        pdf.close()  
        matplotlib.pyplot.close("all")  
         
    return moments 
           
          
