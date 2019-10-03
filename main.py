#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 19:14:08 2019

@author: Egor Kozlov
"""

import numpy as np
import matplotlib.pyplot as plt
from timeit import default_timer
from setup import ModelSetup


#from platform import system


#if system() != 'Darwin':
from solver_couples import vm_period_zero_grid_massive as vm_period_zero_grid
#else:
#    from solver_couples import vm_period_zero_grid_loop as vm_period_zero_grid



if __name__ == '__main__':
    
    
    
    
    
    plt.figure()
    start = default_timer()
    
    
    setup = ModelSetup()
    
    
    print('Setting up done at {}'.format(default_timer()-start))
    
    
    # compute terminal values
    
    Vval_postren, VFval_postren, VMval_postren = setup.vm_last_grid()
    VFval_single, VMval_single = setup.vs_last_grid()
    
    print('Computing on grid done at {}'.format(default_timer()-start))
    
    # assemble them to V structure
    V_last = { 'M':  {'V':Vval_postren,'VF':VFval_postren,'VM':VMval_postren},
          'SF': {'V':VFval_single},
          'SM': {'V':VMval_single}
         }
    
    
    from integrator_singles import ev_after_savings_grid_all_z
    from renegotiation import v_last_period_renegotiated
    from integrator_couples import ev_couple_after_savings
    
    
    
    solution = list()    
    V = V_last
    
    from solver_singles import v_period_zero_grid
    
    for female in [True,False]:
        
        EV_integrated = ev_after_savings_grid_all_z(setup,V,setup.agrid,female)
        print('Integration done at {}'.format(default_timer()-start))        
        
        solution.append(v_period_zero_grid(setup,setup.agrid,EV_integrated,female))
        print('Optimization singles for period 0 done at {}'.format(default_timer()-start))
        
    evc = ev_couple_after_savings(setup,V_last)
    evc = tuple(np.float32(x) for x in evc) # type conversion
    
    print('Integration for couples done at {}'.format(default_timer()-start))
    
    (V_couple,c_couple,s_couple) = vm_period_zero_grid(setup,setup.agrid,evc)
    print('Optimization for couples done at {}'.format(default_timer()-start))
    
    # one more period
    
    
    V_more = {'M': V_couple, 'SF': {'V':solution[0][0]}, 'SM': {'V':solution[1][0]}}
    
    print('Renegotiation for couples done at {}'.format(default_timer()-start))
    
    
    solution_more = list()
    for female in [True,False]:
        
        EV_integrated = ev_after_savings_grid_all_z(setup,V_more,setup.agrid,female)
        print('Integration done at {}'.format(default_timer()-start))        
        
        solution_more.append(v_period_zero_grid(setup,setup.agrid,EV_integrated,female))
        print('Optimization singles for period 0 done at {}'.format(default_timer()-start))
    
    
    evc_more = ev_couple_after_savings(setup,V_more)    
    evc_more = tuple(np.float32(x) for x in evc_more) # type conversion
    
    print('Integration for couples done at {}'.format(default_timer()-start))
    
    
    (V_couple_more,c_couple_more,s_couple_more) = vm_period_zero_grid(setup,setup.agrid,evc_more)
    print('Optimization for couples done at {}'.format(default_timer()-start))
    
    
    
    plt.plot(setup.agrid,solution_more[0][2][:,4],label="female, 0")
    plt.plot(setup.agrid,solution_more[1][2][:,4],label="male, 0")
    
    plt.plot(setup.agrid,solution[0][2][:,4],label="female, 1")
    plt.plot(setup.agrid,solution[1][2][:,4],label="male, 1")
    
    
    plt.legend()
    print('Time elapsed is {}'.format(default_timer()-start))
    

