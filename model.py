#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This collects routines for setting up the model.

Created on Tue Oct 29 17:01:07 2019

@author: egorkozlov
"""


#from platform import system

import numpy as np
from timeit import default_timer




#if system() != 'Darwin':
from setup import ModelSetup
from simulations import Agents
from solver_couples import vm_period_zero_grid_massive as vm_period_zero_grid
from solver_singles import v_period_zero_grid
from integrator_singles import ev_after_savings_grid_all_z
from integrator_couples import ev_couple_after_savings


class Model(object):
    def __init__(self,iterator_name='default',**kwargs):
        self.setup = ModelSetup(**kwargs)
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
        
    def time_statistics(self,remove_worst=True,remove_single=False):
        for what, timelist in self.time_dict.items():
            
            if remove_single and len(timelist) == 1: continue
            
            time_arr = np.array(timelist)
            
            extra = ''
            if remove_worst and time_arr.size > 1:
                time_worst = time_arr.max()
                time_arr = time_arr[time_arr<time_worst]
                extra = ' (excl the worst)'
                
            av_time = round(np.mean(time_arr),2)            
            print('On average {} took {} sec{}'.format(what,av_time,extra))
            
    
    def _get_iterator(self,name='default'):
        # this thing returns two functions: iterate and initialize
        # it can put a timer inside of them
        # it can also do many other things potentially
        
        
        # this is not the best organization but kind of works
        # this allows to use different methods for iteration/initialization
        # as long as they all are specified here or imported
        
        # first we define the iterator
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
            
        # and the integrator   
        def v_integrator(setup,desc,t,V_next):
            
            if desc == 'Female, single' or desc == 'Male, single':
                female = (desc == 'Female, single')
                EV = ev_after_savings_grid_all_z(setup,V_next,setup.agrid,female,t)
            elif desc == 'Couple':
                EV = ev_couple_after_savings(setup,V_next,t)
            return EV
            
        
        
        # then we wrap them into two routines  
        
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
            
            for desc in self.setup.state_names:
                if t == T-1:
                    V_d = self.initializer(desc,t)
                else:
                    V_d = self.iterator(desc,t,Vnext)                    
                Vnow.update(V_d)
            
            self.V = [Vnow] + self.V
            
            
            
            
    def solve_sim(self):
        self.solve()
        self.agents = Agents(self)
        self.agents.simulate()
        
            
    def solve_marriage(self,sf,sm,izf,izm,ipsi,t=0):
        # this is a legacy code that is meant for comparison of unilateral
        # and bilateral divorce
        
        # sf, sm, izf, izm are 1-dim np.arrays 
        # the result is of shape (sf.size,izf.size) (so it is computed for all
        # combinations of sf/sm and izf/izm/ipsi)
        
        # so far it is stable only when izf/izm/ipsi and sf and sm have more 
        # than one element due to shape issues
        
        # this is value of a potential couple that is about to enter period t
        
        
        inds_tuple = (izf,izm,ipsi)
        V = self.V[t]
        Vout_f, Vout_m, ismar, thetaout, technical = \
            v_mar(self.setup,V,sf,sm,inds_tuple,return_all=True,combine=False)
        return Vout_f, Vout_m, ismar, thetaout
    