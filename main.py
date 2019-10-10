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
from solver_singles import v_period_zero_grid
from integrator_singles import ev_after_savings_grid_all_z
from integrator_couples import ev_couple_after_savings
#else:
#    from solver_couples import vm_period_zero_grid_loop as vm_period_zero_grid


class Model(object):
    def __init__(self,iterator_name='default'):
        self.setup = ModelSetup()
        self.iterator, self.initializer = self._get_iterator(iterator_name)
        self.start = default_timer()
        self.last = default_timer()
        self.time_dict = dict()
        
    def time(self,whatisdone):
        
        total_time = default_timer() - self.start
        last_time = default_timer() - self.last
        
        def r(x): return round(x,2)
        
        print('{} is done in {} sec, total {} sec'.format(whatisdone,r(last_time),r(total_time)))
        self.last = default_timer()
    
        if whatisdone in self.time_dict:
            self.time_dict[whatisdone] = self.time_dict[whatisdone] + [last_time]
        else:
            self.time_dict[whatisdone] = [last_time]
        
    def time_statistics(self,remove_worst=True,remove_single=True):
        for what, timelist in self.time_dict.items():
            
            if remove_single and len(timelist) == 1: continue
            
            time_arr = np.array(timelist)
            
            extra = ''
            if remove_worst:
                time_worst = time_arr.max()
                time_arr = time_arr[time_arr<time_worst]
                extra = ' (excl the worst)'
                
            av_time = round(np.mean(time_arr),2)            
            print('On average {} took {} sec{}'.format(what,av_time,extra))
            
            
            
    
    
    
    def _get_iterator(self,name='default'):
        
        if name == 'default' or name == 'default-timed':
            timed = (name == 'default-timed')
            def iterate(desc,t,Vnext):
                EV = v_integrator(self.setup,desc,t,Vnext)
                if timed: self.time('Integration for {}'.format(desc))
                vout = v_iterator(self.setup,desc,t,EV)
                if timed: self.time('Optimization for {}'.format(desc))
                return vout
            def initialize(desc,t):
                vout = v_iterator(self.setup,desc,None)
                if timed: self.time('Initialization for {}'.format(desc))
                return vout
        else:
            raise Exception('unsupported name')
            
        return iterate, initialize
        
        
    def solve(self):
        T = self.setup.pars['T']
        self.V = list()
        
        
        for t in reversed(range(T)):
            Vnow = dict()
            
            Vnext = self.V[0] if t<T-1 else None
            
            for desc in ['Female, single','Male, single','Couple']:
                if t == T-1:
                    V_d = self.initializer(desc,t)
                else:
                    V_d = self.iterator(desc,t,Vnext)                    
                Vnow.update(V_d)
            
            self.V = [Vnow] + self.V
            
            
        
# this is not the best organization but kind of works
def v_iterator(setup,desc,t,EV=None):
    # this takes integrated future type-specific value function and returns
    # this period value function. Integration is done separately.
    # If None is feeded for EV this assumes that we are in the last period
    # and returns last period value
    
    if desc == 'Female, single' or desc == 'Male, single':
        female = (desc == 'Female, single')
        if EV is None:            
            V, c, s = setup.vs_last_grid(female,return_cs=True)
        else:
            V, c, s = v_period_zero_grid(setup,setup.agrid,EV,female)             
        return {desc: {'V':V,'c':c,'s':s}}   
     
    elif desc == 'Couple':
        if EV is None:
            V, VF, VM, c, s = setup.vm_last_grid(return_cs=True)
        else:
            V, VF, VM, c, s = vm_period_zero_grid(setup,setup.agrid,EV)            
        return {desc: {'V':V,'VF':VF,'VM':VM,'c':c,'s':s}}
    
    
def v_integrator(setup,desc,t,V_next):
    
    if desc == 'Female, single' or desc == 'Male, single':
        female = (desc == 'Female, single')
        EV = ev_after_savings_grid_all_z(setup,V_next,setup.agrid,female,t)
    elif desc == 'Couple':
        EV = ev_couple_after_savings(setup,V_next,t)
    return EV
    
        


if __name__ == '__main__':
    
    mdl = Model(iterator_name='default-timed')
    mdl.solve()
    mdl.time_statistics()
    
    

