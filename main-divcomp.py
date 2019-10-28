#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 19:14:08 2019

@author: Egor Kozlov
"""

import numpy as np
from main import Model
        


if __name__ == '__main__':
    
    mdl = Model(iterator_name='default-timed',
                divorce_costs={'unilateral_divorce':True},T=3)
    mdl.solve()
    mdl.time_statistics()
    
    s = mdl.setup
    
    
    sf = np.array([4.0,8.0])
    sm = np.array([8.0,4.0])
    
    t = 0 # at what
    
    izf = np.array( [4,5], dtype=np.int32) # z position out of setup.pars['nf']
    izm = np.array( [3,4], dtype=np.int32) # z position
    ipsi = np.array([7,12],dtype=np.int32) # psi
    
    print( 'zf value is {}'.format( s.exogrid.zf_t[t][izf]) )
    print( 'zm value is {}'.format( s.exogrid.zm_t[t][izm]) )
    print('psi value is {}'.format( s.exogrid.psi_t[t][ipsi]) )
    
    Vout_f, Vout_m, ismar, thetaout = mdl.solve_marriage(sf,sm,izf,izm,ipsi)
    print(thetaout)
    
    
    mdl2 = Model(iterator_name='default-timed',
                 divorce_costs={'unilateral_divorce':False},T=3)
    
    mdl2.solve()
    Vout_f2, Vout_m2, ismar2, thetaout2 = mdl2.solve_marriage(sf,sm,izf,izm,ipsi)
    print(thetaout2)
