#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 14 10:58:43 2019

@author: Egor
"""


# this defines model residuals


import numpy as np
xdef = np.array([0.05,0.01,0.02,0.7,0.25])


# return format is any combination of 'distance', 'all_residuals' and 'model'
# we can add more things too for convenience
def mdl_resid(x=xdef,return_format=['distance'],verbose=False,calibration_report=False):
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
                  verbose_sim=verbose)
    
    Tret = mdl.setup.pars['Tret']
    
    
    
    
    mean_mar = np.mean(mdl.moments['share mar'][1:Tret])
    mean_coh = np.mean(mdl.moments['share coh'][1:Tret])
    
    marcoh_ratio = mean_mar / mean_coh
    
    fls_ratio = np.mean(mdl.moments['flsm'][1:Tret])/np.mean(mdl.moments['flsc'][1:Tret])
    
    haz_sep = mdl.moments['hazard sep'][0]
    haz_div = mdl.moments['hazard div'][0]
    haz_mar = mdl.moments['hazard mar'][0]
    
    coh_ret = mdl.moments['share coh'][Tret-1]
    mar_ret = mdl.moments['share mar'][Tret-1]
    
    resid = [0.0]
    resid += [(coh_ret - 0.1)**2]
    resid += [(mar_ret - 0.8)**2]
    resid += [((fls_ratio - 0.8)**2)]
    resid += [((haz_mar - 0.1)**2)]
    resid += [(haz_sep - 0.15)**2]
    resid += [(haz_div - 0.05)**2]
    
    def distance_pso(particle):
        return 
    
    resid = [v if not (np.isnan(v) or np.isinf(v)) else 1.0e6 for v in resid]
    out = np.sum(resid)
    
    
    if calibration_report:
        print('')
        print('')
        print('Calibration report')
        print('ulost = {:.4f} , s_psi = {:.4f}, s_psi0 = {:.4f}, uls = {:.4f}, pmeet = {:.4f}'.format(ulost,sigma_psi,sigma_psi_init,uls, pmeet))
        print('')
        print('')
        print('At retirement {:.4f} mar and {:.4f} cohab'.format(mar_ret,coh_ret))
        print('All-t ratio of marriages to cohabitation is {:.4f}'.format(marcoh_ratio))
        print('Hazard of sep is {:.4f}, hazard of div is {:.4f}'.format(haz_sep,haz_div))        
        print('Hazard of Marriage is {:.4f}'.format(haz_mar))
        print('Calibration residual is {:.4f}'.format(out))
        print('')
        print('')
        print('End of calibration report')
        print('')
        print('')
    
    
    
    out_dict = {'distance':out,'all residuals':resid,'model':mdl}
    out = [out_dict[key] for key in return_format]
    return out
