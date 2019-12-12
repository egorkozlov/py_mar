#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
aCreated on Tue Sep 17 19:14:08 2019

@author: Egor Kozlov
"""





if __name__ == '__main__':
    
    try:
        from IPython import get_ipython
        get_ipython().magic('reset -f')
    except:
        pass


from platform import system
    
if system() != 'Darwin' and system() != 'Windows':   
    import os
    os.environ['QT_QPA_PLATFORM']='offscreen'


import numpy as np
from numpy.random import random_sample as rs
from tiktak import tiktak

xdef = np.array([0.05,0.01,0.02,0.7,0.25])
def mdl_resid(x=xdef):
    from model import Model
    from setup import DivorceCosts
    
    
    ulost = x[0]
    sigma_psi = x[1]
    sigma_psi_init = x[2]
    pmeet = x[3]#np.exp(x[3])/(1+np.exp(x[3]))
    uls = x[4]
    

    dc = DivorceCosts(unilateral_divorce=True,assets_kept = 1.0,u_lost_m=ulost,u_lost_f=ulost,eq_split=0.0)
    sc = DivorceCosts(unilateral_divorce=True,assets_kept = 1.0,u_lost_m=0.00,u_lost_f=0.00)
    mdl = Model(iterator_name='default',divorce_costs=dc,
                separation_costs=sc,sigma_psi=sigma_psi,
                sigma_psi_init=sigma_psi_init,
                pmeet=pmeet,uls=uls)
    
    mdl.solve_sim(simulate=True)
    
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
    resid += [((marcoh_ratio - 1.1)**2)*(marcoh_ratio<1.1)]
    resid += [((fls_ratio - 0.8)**2)*(fls_ratio > 0.8)]
    resid += [((haz_mar - 0.15)**2)]
    resid += [(haz_sep - 0.2)**2]
    resid += [(haz_div - 0.05)**2]
    
    def distance_pso(particle):
        return 
    
    resid = [v if not (np.isnan(v) or np.isinf(v)) else 1.0e6 for v in resid]
    out = np.sum(resid)
    
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
    
    
    
    
    return out



if __name__ == '__main__':
    
    
    #If on server set Display
    
 
   
            
    #Create grids of parameters
    sigma_psi_g=np.linspace(0.01,0.3,3)
    sigma_psi_init_g=np.linspace(0.05,0.5,3)
    di_co_g=np.linspace(0.05,0.3,3)
    bila=np.array([False,True])
    
    
    
    #Initialize the file with parameters


    x0 = np.array([0.05,0.01,0.02,0.7,0.25])
    lb= np.array([0.0,0.005,0.015,0.4,0.01])
    ub= np.array([1.0,1.0,0.45,1.0,0.4])
    ub[3]=min(ub[3],1.0)
    
    
    ##### FIRST LET'S TRY TO RUN THE FUNCTION IN FEW POINTS
    
    print('Testing the workers...')
    from p_client import compute_for_values
    pts = [lb + rs(lb.shape)*(ub-lb) for _ in range(2)]
    pts = [('compute',x) for x in pts]    
    outs = compute_for_values(pts,timeout=3600.0)
    print('Everything worked, output is {}'.format(outs))
    
    
    
    x0 = np.array([0.05,0.01,0.02,0.7,0.25])
    lb= x0*0.5#np.array([0.01,0.005,0.015,0.4,0.01])
    ub= x0*2.0#np.array([0.1,0.3,0.45,1.0,0.4])
    ub[3]=min(ub[3],1.0)
   
    print('')
    print('')
    print('running tic tac...')
    print('')
    print('')
    
   

    #Tik Tak Optimization
    param=tiktak(1,2,1,lb,ub,mdl_resid,tole=1e-3,nelder=False,refine=False)
    
    print('f is {} and x is {}'.format(param[0],param[1]))
    
  
    
    
    xlist = [lb + rs(x0.shape)*(ub-lb) for _ in range(50)]
    from p_client import compute_for_values
    out = compute_for_values(xlist)
    
    ibest = np.argmin(out)
    xbest = xlist[ibest]
    print(out[ibest],xbest)
    
   
        

