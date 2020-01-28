#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 14 10:58:43 2019

@author: Egor
"""


# this defines model residuals
import numpy as np
import pickle, dill
import os

xdef = np.array([0.05,0.01,0.02,0.7,0.25,0.0001,0.5,0.5])


# return format is any combination of 'distance', 'all_residuals' and 'models'
# we can add more things too for convenience
def mdl_resid(x=xdef,save_to=None,load_from=None,return_format=['distance'],
              solve_transition=False,
              store_path = None,
              verbose=False,calibration_report=False,draw=False,graphs=False,
              rel_diff=True):
    
    
    
    from model import Model
    from setup import DivorceCosts
    from simulations import Agents
    from moments import moment
 
    ulost = x[0] #min(x[0],1.0)
    mshift=x[5]
    sigma_psi = x[1] # max(x[1],0.00001)
    sigma_psi_init = x[1]*x[2] # max(x[2],0.00001) # treat x[2] as factor
    pmeet = x[3] # #min(x[3],1.0)#np.exp(x[3])/(1+np.exp(x[3]))
    pls = x[6] #max(min(x[6],1.0),0.0)
    util_alp = x[4]
    util_kap = x[7]
    
    
    
    # this is for the default model
    dc = DivorceCosts(unilateral_divorce=True,assets_kept = 1.0,u_lost_m=ulost,u_lost_f=ulost,eq_split=0.0)
    sc = DivorceCosts(unilateral_divorce=True,assets_kept = 1.0,u_lost_m=0.00,u_lost_f=0.00)
    
    
    
    
    iter_name = 'default' if not verbose else 'default-timed'
    
    
    def join_path(name,path):
        return os.path.join(path,name)
    
    
    
    
    
    
    if load_from is not None:
        if type(load_from) is not list:
            load_from = [load_from]
        if store_path is not None:
            load_from = [join_path(n,store_path) for n in load_from]
    
    
    if save_to is not None:
        if type(save_to) is not list:
            save_to = [save_to]
        if store_path is not None:
            save_to = [join_path(n,store_path) for n in save_to]
    
    
                
    
    
    if load_from is None:
        
        kwords = dict(sigma_psi=sigma_psi,
                        sigma_psi_init=sigma_psi_init,
                        pmeet=pmeet,util_alp=util_alp,util_kap=util_kap,
                        pls=pls,u_shift_mar=mshift)
        
        if not solve_transition:
            
            mdl = Model(iterator_name=iter_name,divorce_costs=dc,
                        separation_costs=sc,**kwords)
            mdl_list = [mdl]
            
        else:
            # specify the changes here manually
            dc_before = DivorceCosts(unilateral_divorce=False,assets_kept = 1.0,u_lost_m=ulost,u_lost_f=ulost,eq_split=0.0)
            dc_after  = DivorceCosts(unilateral_divorce=True,assets_kept = 1.0,u_lost_m=ulost,u_lost_f=ulost,eq_split=0.0)
            
            
            
            
            mdl_before = Model(iterator_name=iter_name,divorce_costs=dc_before,
                        separation_costs=sc,**kwords)
            
            mdl_after = Model(iterator_name=iter_name,divorce_costs=dc_after,
                        separation_costs=sc,**kwords)  
            
            mdl = mdl_after # !!! check if this makes a difference
            # I think that it is not used for anything other than getting 
            # setup for plotting
            
            mdl_list = [mdl_before,mdl_after]
            
    else:       
        mdl_list = [dill.load(open(l,'rb+')) for l in load_from]
        mdl = mdl_list[0]
        
        if solve_transition:
            if len(mdl_list) < 2:
                print('Warning: you supplied only one model, so no transition is computed')
    
    if save_to is not None:
        
        if not solve_transition:
            if len(save_to) > 1:
                print('warning: too much stuff is save_to')
            dill.dump(mdl,open(save_to[0],'wb+'))            
            
        else:            
            if len(save_to) > 1:
                [dill.dump(m_i,open(n_i,'wb+')) 
                    for (m_i,n_i) in zip(mdl_list,save_to)]
            else:
                print('Warning: file names have change to write two models, \
                      please provide the list of names if you do not want this')
                dill.dump(mdl_before,open(save_to[0] + '_before','wb+'))
                dill.dump(mdl_after, open(save_to[0] + '_after','wb+'))
                
    
    ##############################################################
    # Build Markov transition processes for models from the data
    #############################################################
    
    #Import Data
    with open('age_uni.pkl', 'rb') as file:
        age_uni=pickle.load(file)
    
    #Transfrom distribution of age at Unilateral Divorce into conditional Probabilities
    #The Probability of Transitioning from Unilateral to Bilateral is always zero
    
    #Transformation of age at uni from actual age to model periods
    change=-np.ones(1000,np.int32)#the bigger is the size of this array, the more precise the final distribution
   
    summa=0.0
    summa1=0.0
    for i in age_uni:
        summa+=age_uni[i]
        change[int(summa1*len(change[:])/sum(age_uni.values())):int(summa*len(change[:])/sum(age_uni.values()))]=(i-18)/mdl.setup.pars['py']
        summa1+=age_uni[i]
    change=np.sort(change, axis=0) 
    
    #Now we compute the actual conditional probabilities
    transition_matrices=list()
    
    #First period treated differently
    pr=np.sum(change<=0)/(np.sum(change<=np.inf))
    transition_matrices=transition_matrices+[np.array([[1-pr,pr],[0,1]])]
    for t in range(mdl.setup.pars['T']-1):
        pr=np.sum(change==t+1)/(np.sum(change<=np.inf))
        transition_matrices=transition_matrices+[np.array([[1-pr,pr],[0,1]])]
        
    
    agents = Agents( mdl_list ,pswitchlist=transition_matrices,verbose=verbose)
    moments = moment(mdl,agents,draw=draw)
    
    ############################################################
    #Build data moments and compare them with simulated ones
    ###########################################################
    
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
    W=packed_data['W']
    dat=np.concatenate((hazm_d,hazs_d,hazd_d,mar_d,coh_d,fls_d,beta_unid_d),axis=0)
    

    #Get Simulated Data
    Tret = mdl.setup.pars['Tret']
    hazm_s = moments['hazard mar'][0:len(hazm_d)]
    hazs_s = moments['hazard sep'][0:len(hazs_d)]
    hazd_s = moments['hazard div'][0:len(hazd_d)]
    mar_s = moments['share mar'][0:len(mar_d)]
    coh_s = moments['share coh'][0:len(coh_d)]
    beta_unid_s=np.ones(1)*moments['beta unid']
    fls_s = np.ones(1)*np.mean(moments['flsm'][1:Tret])/np.mean(moments['flsc'][1:Tret])
    sim=np.concatenate((hazm_s,hazs_s,hazd_s,mar_s,coh_s,fls_s,beta_unid_s),axis=0)



    if len(dat) != len(sim):
        sim = np.full_like(dat,1.0e6)
        
   
    res_all=(dat-sim)
  
    
    
    if verbose:
        print('data moments are {}'.format(dat))
        print('simulated moments are {}'.format(sim))
    
    resid_all = np.array([x if (not np.isnan(x) and not np.isinf(x)) else 1e6 for x in res_all])
    
    
    
    resid_sc = resid_all*np.sqrt(np.diag(W)) # all residuals scaled
    
    dist = np.dot(np.dot(resid_all,W),resid_all)


    print('Distance is {}'.format(dist))
    
    
    
    if calibration_report:
        print('')
        print('')
        print('Calibration report')
        print('ulost = {:.4f} , s_psi = {:.4f}, s_psi0 = {:.4f}, uls = {:.4f}, pmeet = {:.4f}'.format(ulost,sigma_psi,sigma_psi_init,uls, pmeet))
        print('')
        print('')
        print('Average {:.4f} mar and {:.4f} cohab'.format(np.mean(mar_s),np.mean(coh_s)))
        print('Hazard of sep is {:.4f}, hazard of div is {:.4f}'.format(np.mean(hazs_s),np.mean(hazd_s)))        
        print('Hazard of Marriage is {:.4f}'.format(np.mean(hazm_s)))
        print('Calibration residual is {:.4f}'.format(dist))
        print('')
        print('')
        print('End of calibration report')
        print('')
        print('')
    
    
    
    out_dict = {'distance':dist,'all residuals':resid_all,
                'scaled residuals':resid_sc,'models':mdl_list,'agents':agents}
    out = [out_dict[key] for key in return_format]
    
  
            
    return out
