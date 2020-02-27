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




# return format is any combination of 'distance', 'all_residuals' and 'models'
# we can add more things too for convenience
def mdl_resid(x=None,save_to=None,load_from=None,return_format=['distance'],
              solve_transition=False,
              store_path = None,
              verbose=False,calibration_report=False,draw=False,graphs=False,
              rel_diff=True):
    
    
    
    from model import Model
    from setup import DivorceCosts
    from simulations import Agents, AgentsPooled
    from moments import moment
    
    
    
    if type(x) is dict:
        params = x
    else:
        from calibration_params import calibration_params     
        lb, ub, x0, keys, translator = calibration_params()        
        if verbose: print('Calibration adjusts {}'.format(keys))        
        if x is None: x = x0
        params = translator(x) # this converts x to setup parameters according to what is in calibration_params
    
    
    
    
    
    
    try:
        ulost = params.pop('ulost') # ulost does not belong to setup parameters so this pop removes it
    except:
        ulost = 0.0
        
    try:
        alost = params.pop('alost')
    except:
        alost = 0.0
    
    # this is for the default model
    dc = DivorceCosts(unilateral_divorce=False,assets_kept = 0.8,u_lost_m=ulost,u_lost_f=ulost,eq_split=1.0)
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
        
        
        
        if not solve_transition:
            
            mdl = Model(iterator_name=iter_name,divorce_costs=dc,
                        separation_costs=sc,**params)
            mdl_list = [mdl]
            
        else:
            # specify the changes here manually
            dc_before = DivorceCosts(unilateral_divorce=False,assets_kept = 0.8,u_lost_m=ulost,u_lost_f=ulost,eq_split=1.0)
            dc_after  = DivorceCosts(unilateral_divorce=True,assets_kept = 0.8,u_lost_m=ulost,u_lost_f=ulost,eq_split=1.0)
            
            
            
            mdl_before = Model(iterator_name=iter_name,divorce_costs=dc_before,
                        separation_costs=sc,**params)
            
            mdl_after = Model(iterator_name=iter_name,divorce_costs=dc_after,
                        separation_costs=sc,**params)  
            
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
        print(i)
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
        
    
    agents_fem = Agents( mdl_list ,female=True,pswitchlist=transition_matrices,verbose=verbose)
    agents_mal = Agents( mdl_list ,female=False,pswitchlist=transition_matrices,verbose=verbose)
    agents_pooled = AgentsPooled([agents_fem,agents_mal])
    
    
    moments = moment(mdl,agents_fem,draw=draw)
    
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
    mean_fls_d=np.ones(1)*packed_data['mean_fls']
    W=packed_data['W']
    dat=np.concatenate((hazm_d,hazs_d,hazd_d,mar_d,coh_d,fls_d,beta_unid_d,mean_fls_d),axis=0)
    

    #Get Simulated Data
    Tret = mdl.setup.pars['Tret']
    hazm_s = moments['hazard mar'][0:len(hazm_d)]
    hazs_s = moments['hazard sep'][0:len(hazs_d)]
    hazd_s = moments['hazard div'][0:len(hazd_d)]
    mar_s = moments['share mar'][0:len(mar_d)]
    coh_s = moments['share coh'][0:len(coh_d)]
    beta_unid_s=np.ones(1)*moments['beta unid']
    mean_fls_s=np.ones(1)*moments['mean_fls']
    fls_s = np.ones(1)*moments['fls_ratio']
    sim=np.concatenate((hazm_s,hazs_s,hazd_s,mar_s,coh_s,fls_s,beta_unid_s,mean_fls_s),axis=0)



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
    
    
    
    del(out_dict)
    
    # memory management
    if 'models' not in return_format:
        for m in mdl_list:
            del(m)
        del mdl_list
        
    if 'agents' not in return_format:
        del(agents)
        
  
            
    return out
