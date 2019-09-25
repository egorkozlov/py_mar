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




if __name__ == '__main__':
    
    
    run_exact = False
    run_exact_with_int = False
    
    
    plt.figure()
    start = default_timer()
    
    
    setup = ModelSetup()
    
    # relevant for plotting
    a0 = 2
    zgrid = setup.exogrid.zf_t[0]
    z0 = setup.exogrid.zf_t[0][4]
    
    print('Setting up done at {}'.format(default_timer()-start))
    
    if run_exact:
        from main_exact import exact_solution
        savings_rate = exact_solution(setup,setup.agrid,np.array([z0]),ireturn=2)
        plt.plot(setup.agrid,savings_rate,label="exact") # these are exact savings on the grid
        print('Exact solution done at {}'.format(default_timer()-start))
    
    
    if run_exact_with_int:
        from main_exact import exact_solution_with_interpolated_ev
        savings_rate_int = exact_solution_with_interpolated_ev(setup,setup.agrid,np.array([z0]),ireturn=2)
        plt.plot(setup.agrid,savings_rate_int,label="interpolation-dumb")
        print('Exact solution with interpolated EV done at {}'.format(default_timer()-start))
    
    
    # compute terminal values
    Vval_postren, VFval_postren, VMval_postren = setup.vm_last_grid()
    VFval_single, VMval_single = setup.vs_last_grid()
    
    print('Computing on grid done at {}'.format(default_timer()-start))
    
    # assemble them to V structure
    V0 = { 'M':  {'V':Vval_postren,'VF':VFval_postren,'VM':VMval_postren},
          'SF': {'V':VFval_single},
          'SM': {'V':VMval_single}
         }
    
    from integrator_singles import ev_after_savings_grid_all_z
    from renegotiation import v_last_period_renegotiated
    from integrator_couples import ev_couple_after_savings
    
    vv = v_last_period_renegotiated(setup,V0)
    print('Renegotiation for period 1 done at {}'.format(default_timer()-start))
    
    
    V1 = { 'M':  {'V':vv[0],'VF':vv[1],'VM':vv[2]},
          'SF': {'V':VFval_single},
          'SM': {'V':VMval_single}
         }
    
    
    evc = ev_couple_after_savings(setup,V1['M'])
    print('Integration for couples done at {}'.format(default_timer()-start))
    
    female = True
    
    solution = list()
    
    V = V0
    
    for female in [True,False]:
        
        EV_integrated = ev_after_savings_grid_all_z(setup,V,setup.agrid,female)
        print('Integration done at {}'.format(default_timer()-start))
        
        from solver_singles import v_period_zero_grid
        solution.append(v_period_zero_grid(setup,setup.agrid,EV_integrated,female))
        print('Optimization singles for period 0 done at {}'.format(default_timer()-start))
        
    
    plt.plot(setup.agrid,solution[0][2][:,4],label="interpolation-smart-0")
    plt.plot(setup.agrid,solution[1][2][:,4],label="interpolation-smart-1")
    
    
    plt.legend()
    print('Time elapsed is {}'.format(default_timer()-start))





