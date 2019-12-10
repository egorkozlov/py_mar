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
#from numba import njit, vectorize
#from memory_profiler import profile
#from IPython import get_ipython
from asizeof import asizeof
import os
import psutil


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
        self.mstart = self.get_mem()
        self.mlast = self.get_mem()
        
        self.setup = ModelSetup(**kwargs)
        self.iterator, self.initializer = self._get_iterator(iterator_name)
        self.start = default_timer()
        self.last = default_timer()
        
        self.time_dict = dict()
        
    def get_mem(self):
        return psutil.Process(os.getpid()).memory_info().rss/1e6
        
        
    def time(self,whatisdone):
        
        total_time = default_timer() - self.start
        last_time = default_timer() - self.last
        
        total_mem = self.get_mem()
        last_mem  = self.get_mem() - self.mlast
        
        def r(x): return round(x,2)
        
        print('{} is done in {} sec, total {} sec, memory used is {} Mb'.format(whatisdone,r(last_time),r(total_time),r(total_mem)))
        self.last = default_timer()
        self.mlast = self.get_mem()
        
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
            #get_ipython().magic('reset -sf')
            if desc == 'Female, single' or desc == 'Male, single':
                female = (desc == 'Female, single')
                if EV is None:            
                    V, c, s = setup.vs_last_grid(female,return_cs=True)
                else:
                    V, c, s = v_iter_single(setup,t,EV,female)             
                return {desc: {'V':V,'c':c,'s':s}}   
             
            elif desc== 'Couple, M' or desc == 'Couple, C':
                if EV is None:
                    V, VF, VM, c, s, fls, V_all_l = setup.vm_last_grid()
                else:
                    V, VF, VM, c, s, fls, V_all_l = v_iter_couple(setup,t,EV)            
                return {desc: {'V':V,'VF':VF,'VM':VM,'c':c,'s':s,'fls':fls,'V_all_l':V_all_l}}
          
            
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
                
                self.wrap_decisions(desc,dec,vout)
                
                return vout, dec
            def initialize(desc,t):
                vout = v_iterator(self.setup,desc,None)
                if timed: self.time('Initialization for {}'.format(desc))
                dec = {}
                self.wrap_decisions(desc,dec,vout)
                return vout, dec
        else:
            raise Exception('unsupported name')
            
            
            
        return iterate, initialize
    
    def wrap_decisions(self,desc,dec,vout):
        # This interpolates consumption, savings and labor supply decisions
        # on fine grid for theta that is used for integration and simulations.
        
        v = vout[desc]
        if desc == 'Couple, M' or desc == 'Couple, C':
            
            #cint = self.setup.v_thetagrid_fine.apply(v['c'],axis=2)
            sint = self.setup.v_thetagrid_fine.apply(v['s'],axis=2).astype(np.float32)
            
            Vint = self.setup.v_thetagrid_fine.apply(v['V_all_l'],axis=2).astype(np.float32)
            
            fls = self.setup.v_thetagrid_fine.apply(Vint,axis=2).argmax(axis=3).astype(np.int8)
            
            dec.update({'s':sint,'fls':fls})
            del sint,fls
        else:
            dec.update({'s':v['s']})
            del v
        
    
    def solve(self,show_mem=False):
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
            

            
            if show_mem:
                print('The size of V is {} giga'.format(asizeof(self.V)/1000000000))
                print('The size of decisions is {} giga'.format(asizeof(self.decisions)/1000000000))
                
        #del self.V,self.decisions    
           
    def solve_sim(self,simulate=True,show_mem=False,draw_moments=True,verbose_sim=False):

        #Solve the model
        self.solve(show_mem=show_mem)
        if not simulate: return
        self.agents = Agents(self,verbose=verbose_sim)
        self.agents.simulate()        
        self.compute_moments(draw=draw_moments)
        
        
    def graph(self,ai,zfi,zmi,psii,ti,thi):        
        #Draw some graph of Value and Policy Functions
        V=graphs(self,ai,zfi,zmi,psii,ti,thi)        
        return V
      
        
    def compute_moments(self,draw=True):
        moment(self,draw=draw)
    
