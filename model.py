#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This collects routines for setting up the model.

Created on Tue Oct 29 17:01:07 2019

@author: egorkozlov
"""


#from platform import system

import numpy as np
import dill as pickle
from timeit import default_timer
import gzip
from numba import njit, vectorize


#if system() != 'Darwin':
from setup import ModelSetup
from graph import graphs
from moments import moment
from simulations import Agents
from solver_couples import v_iter_couple
from solver_singles import v_iter_single
from integrator_singles import ev_single
from integrator_couples import ev_couple_m_c


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
                    V, c, s = v_iter_single(setup,EV,female)             
                return {desc: {'V':V,'c':c,'s':s}}   
             
            elif desc== 'Couple, M' or desc == 'Couple, C':

                if EV is None:
                    V, VF, VM, c, s = setup.vm_last_grid(return_cs=True)
                else:
                    V, VF, VM, c, s = v_iter_couple(setup,EV)            
                return {desc: {'V':V,'VF':VF,'VM':VM,'c':c,'s':s}}
          
            
        # and the integrator   
        def v_integrator(setup,desc,t,V_next):
            
            if desc == 'Female, single' or desc == 'Male, single':
                female = (desc == 'Female, single')
                EV, dec = ev_single(setup,V_next,setup.agrid_s,female,t)
            elif desc == 'Couple, M':
                EV, dec = ev_couple_m_c(setup,V_next,t,True)
            elif desc == 'Couple, C':
                EV, dec = ev_couple_m_c(setup,V_next,t,False)
            return EV, dec
            
        
        
        # then we wrap them into two routines  
        
        if name == 'default' or name == 'default-timed':
            timed = (name == 'default-timed')
            def iterate(desc,t,Vnext):
                EV, dec = v_integrator(self.setup,desc,t,Vnext)
                if timed: self.time('Integration for {}'.format(desc))
                vout = v_iterator(self.setup,desc,t,EV)
                if timed: self.time('Optimization for {}'.format(desc))
                return vout, dec
            def initialize(desc,t):
                vout = v_iterator(self.setup,desc,None)
                if timed: self.time('Initialization for {}'.format(desc))
                return vout, None
        else:
            raise Exception('unsupported name')
            
            
            
        return iterate, initialize
        
        
    def solve(self):
        T = self.setup.pars['T']
        self.V = list()
        self.decisions = list()
        

        
        
        for t in reversed(range(T)):
            Vnow = dict()
            decnow = dict()
            
            Vnext = self.V[0] if t<T-1 else None
            
            for desc in self.setup.state_names:
                if t == T-1:
                    V_d, dec = self.initializer(desc,t)
                else:
                    V_d, dec = self.iterator(desc,t,Vnext)                    
                Vnow.update(V_d)
                decnow.update({desc:dec})
            
            self.V = [Vnow] + self.V
            self.decisions = [decnow] + self.decisions

            
            
            
    def solve_sim(self,simulate=True):

        #Solve the model
        self.solve()
        if not simulate: return
        #Simulate the model
        self.agents = Agents(self)
        self.agents.simulate()
        moment(self.agents,True)
        
        
        #return gassets,iexo,state,gtheta
        
    def graph(self,ai,zfi,zmi,psii,ti,thi):
        
        #Draw some graph of Value and Policy Functions
        V=graphs(self.setup,self.V,self.decisions,ai,zfi,zmi,psii,ti,thi)
        
        return V
      
    