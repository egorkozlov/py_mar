#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 14 10:58:43 2019

@author: Egor
"""


# this defines model residuals
import numpy as np
import pickle
xdef = np.array([0.05,0.01,0.02,0.7,0.25])


# return format is any combination of 'distance', 'all_residuals' and 'model'
# we can add more things too for convenience
def mdl_resid(x=xdef,return_format=['distance'],verbose=False,calibration_report=False,draw=False,graphs=False):
    from model import Model
    from setup import DivorceCosts
    
 
    ulost = x[0]
    sigma_psi = max(x[1],0.00001)
    sigma_psi_init = max(x[2],0.00001)
    pmeet = min(x[3],1.0)#np.exp(x[3])/(1+np.exp(x[3]))
    uls = x[4]
    

    dc = DivorceCosts(unilateral_divorce=True,assets_kept = 1.0,u_lost_m=ulost,u_lost_f=ulost,eq_split=0.0)
    sc = DivorceCosts(unilateral_divorce=True,assets_kept = 1.0,u_lost_m=0.00,u_lost_f=0.00)
    
    
    iter_name = 'default' if not verbose else 'default-timed'
    
    mdl = Model(iterator_name=iter_name,divorce_costs=dc,
                separation_costs=sc,sigma_psi=sigma_psi,
                sigma_psi_init=sigma_psi_init,
                pmeet=pmeet,uls=uls)
    
    mdl.solve_sim(simulate=True,show_mem=verbose,
                  verbose_sim=verbose,draw_moments=draw)
    
    
    ############################################################
    #Build data moments and compare them with simulated ones
    ###########################################################
    
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
    dat=np.concatenate((hazm_d,hazs_d,hazd_d,mar_d,coh_d,fls_d),axis=0)
    W=packed_data[6]

    #Get Simulated Data
    Tret = mdl.setup.pars['Tret']
    hazm_s = mdl.moments['hazard mar'][0:len(hazm_d)]
    hazs_s = mdl.moments['hazard sep'][0:len(hazs_d)]
    hazd_s = mdl.moments['hazard div'][0:len(hazd_d)]
    mar_s = mdl.moments['share mar'][0:len(mar_d)]
    coh_s = mdl.moments['share coh'][0:len(coh_d)]
    fls_s = np.ones(1)*np.mean(mdl.moments['flsm'][1:Tret])/np.mean(mdl.moments['flsc'][1:Tret])
    sim=np.concatenate((hazm_s,hazs_s,hazd_s,mar_s,coh_s,fls_s),axis=0)

    #Get residuals
    vec=dat-sim
    resid = np.dot(np.dot(vec,W),vec)
    
    def distance_pso(particle):
        return 
    
    resid = [resid if not (np.isnan(resid) or np.isinf(resid)) else 1.0e6]
    out = resid
    
    
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
        print('Calibration residual is {:.4f}'.format(out))
        print('')
        print('')
        print('End of calibration report')
        print('')
        print('')
    
    
    
    out_dict = {'distance':out,'all residuals':resid,'model':mdl}
    out = [out_dict[key] for key in return_format]
    
    #For memory reason:delete stuff
    if not draw:
        if not graphs:
            del mdl
            
    return out
