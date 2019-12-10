#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
aCreated on Tue Sep 17 19:14:08 2019

@author: Egor Kozlov
"""


if __name__ == '__main__':
    
    
    #Clean Memory
    
    try:
        from IPython import get_ipython
        get_ipython().magic('reset -f')
    except:
        pass
    
    #If on server set Display
    from platform import system
    
    if system() != 'Darwin' and system() != 'Windows':   
        import os
        os.environ['QT_QPA_PLATFORM']='offscreen'
 
    import numpy as np
    from model import Model
    from setup import DivorceCosts
    from scipy.optimize import minimize
    

            
    #Create grids of parameters
    sigma_psi_g=np.linspace(0.01,0.3,3)
    sigma_psi_init_g=np.linspace(0.05,0.5,3)
    di_co_g=np.linspace(0.05,0.3,3)
    bila=np.array([False,True])
    
    
    
    #Initialize the file with parameters

    import xlwt 
    from xlwt import Workbook 
   
    
    def mdl_resid(x):
        
        ulost = x[0]
        sigma_psi = np.exp(x[1])
        sigma_psi_init = np.exp(x[2])
        

        dc = DivorceCosts(unilateral_divorce=True,assets_kept = 1.0,u_lost_m=ulost,u_lost_f=ulost,eq_split=0.0)
        sc = DivorceCosts(unilateral_divorce=True,assets_kept = 1.0,u_lost_m=0.00,u_lost_f=0.00)
        mdl = Model(iterator_name='default',divorce_costs=dc,separation_costs=sc,sigma_psi=sigma_psi,sigma_psi_init=sigma_psi_init)
        
        mdl.solve_sim(simulate=True)
        
        Tret = mdl.setup.pars['Tret']
        
        
        
        
        mean_mar = np.mean(mdl.moments['share mar'][1:Tret])
        mean_coh = np.mean(mdl.moments['share coh'][1:Tret])
        
        marcoh_ratio = mean_mar / mean_coh
        
        fls_ratio = np.mean(mdl.moments['flsm'][1:Tret])/np.mean(mdl.moments['flsc'][1:Tret])
        
        haz_sep = np.mean(mdl.moments['hazard sep'])
        haz_div = np.mean(mdl.moments['hazard div'])
        
        haz_rat = haz_sep / haz_div
        
        coh_ret = mdl.moments['share coh'][Tret-1]
        mar_ret = mdl.moments['share mar'][Tret-1]
        
        resid = [0.0]
        resid += [(coh_ret - 0.1)**2]
        resid += [(mar_ret - 0.7)**2]
        resid += [((marcoh_ratio - 1.1)**2)*(marcoh_ratio<1.1)]
        resid += [((fls_ratio - 0.8)**2)*(fls_ratio > 0.8)]
        resid += [((haz_rat - 2)**2)*(haz_rat < 2)]
        resid += [(haz_sep - 0.2)**2]
        resid += [(haz_div - 0.05)**2]
        
        print('')
        print('')
        print('Calibration report')
        print('ulost = {} , sigma_psi = {}, sigma_psi_init = {}'.format(ulost,sigma_psi,sigma_psi_init))
        print('')
        print('')
        print('At retirement {} mar and {} cohab'.format(mar_ret,coh_ret))
        print('All-t ratio of marriages to cohabitation is {}'.format(marcoh_ratio))
        print('Hazard of sep is {}, hazard of div is {}'.format(haz_sep,haz_div))        
        print('Relative hazard of sep to div is {}'.format(haz_rat))
        print('')
        print('')
        print('End of calibration report')
        print('')
        print('')
        
        
        
        resid = [v if not (np.isnan(v) or np.isinf(v)) else 1.0e6 for v in resid]
        return np.sum(resid)
    
    
    x0 = np.array([0.1,np.log(0.05),np.log(0.15)])
    print('Initial value is {}'.format(mdl_resid(x0)))
    print('')
    print('')
    print('minimizing...')
    print('')
    print('')
    
    
    res = minimize(mdl_resid,x0,method='Nelder-Mead')
    
    print('x is {} and fun is {}'.format(res.x,res.fun))
    
    print('Final value is {}'.format(mdl_resid(res.x)))
    
    
#    #Indexes for the graphs
#    if graphs:
#        ai=0
#        zfi=0
#        zmi=4
#        psii=5
#        ti=1
#        thi=10
#        
#        #Actual Graphs
#        mdl.graph(ai,zfi,zmi,psii,ti,thi)
#        
#        #If you plan to use graphs only once, deselect below to save space on disk
#        #os.remove('name_model.pkl') 
    


   
        


    