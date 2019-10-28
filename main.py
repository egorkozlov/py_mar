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


from interp_np import interp            
from mc_tools import mc_simulate
from ren_mar import v_mar, v_ren


#from platform import system


#if system() != 'Darwin':
from solver_couples import vm_period_zero_grid_massive as vm_period_zero_grid
from solver_singles import v_period_zero_grid
from integrator_singles import ev_after_savings_grid_all_z
from integrator_couples import ev_couple_after_savings

#else:
#    from solver_couples import vm_period_zero_grid_loop as vm_period_zero_grid


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
    
    


class Agents:
    
    
    def __init__(self,M,N=1000,T=None,verbose=True):
        if T is None:
            T = M.setup.pars['T']
            
            
            
        self.M = M
        self.V = M.V
        self.setup = M.setup
        self.state_names = self.setup.state_names
        self.N = N
        self.T = T
            
        self.verbose = verbose
        
        np.random.seed(18)    
        
        self.timer = M.time
        self._iassets = np.zeros((N,T),np.int32)            
        self._wnassets = np.ones((N,T),np.float32)
        
        self._itheta  = np.zeros((N,T),dtype=np.int32)           
        self._wntheta  = np.zeros((N,T),dtype=np.float32)
        
        
        self.assets = np.zeros((N,T),dtype=np.float32)
        self.theta = np.full((N,T),np.nan,dtype=np.float32)
        self.iexo = np.zeros((N,T),np.int32)
        self.iexo[:,0] = np.random.randint(0,self.setup.pars['n_zf'],size=N)
        #self.iassets[:,0] = np.random.randint(0,100,size=N)
        
        self.state = np.full((N,T),np.nan,dtype=np.int32)
        #self.state[:,0] = np.random.randint(0,2,size=N,dtype=np.int32)
        
        
        
        
        self.state_codes = dict()
        self.has_theta = list()
        for i, name in enumerate(self.setup.state_names):
            self.state_codes[name] = i
            self.has_theta.append((name=='Couple'))
        
        
        self.state[:,0] = 0
        
        self.timer('Simulations, creation')
    
    def simulate(self):
        
        for t in range(self.T-1):
            
            self.anext(t)            
            self.iexonext(t)            
            self.statenext(t)
            self.timer('Simulations, iteration')
        
        
    def anext(self,t):
        
        for ist, sname in enumerate(self.state_codes):
            
            is_state = (self.state[:,t]==ist)
            
            use_theta = self.has_theta[ist]
            
            nst = np.sum(is_state)
            if nst==0:
                continue
            
            ind = np.where(is_state)
            
            agrid = self.setup.agrid
            
            
            if not use_theta:
                s_i  = self.V[t][sname]['s'][self._iassets[ind,t],  self.iexo[ind,t]].reshape(nst)
                s_ip = self.V[t][sname]['s'][self._iassets[ind,t]+1,self.iexo[ind,t]].reshape(nst)
            else:
                wt_p = self._wntheta[ind,t].reshape(nst)                
                # bilinear interpolation                
                s_i_t0 = self.V[t][sname]['s'][self._iassets[ind,t],  self.iexo[ind,t], self._itheta[ind,t]].reshape(nst)
                s_i_tp = self.V[t][sname]['s'][self._iassets[ind,t],  self.iexo[ind,t], self._itheta[ind,t]+1].reshape(nst)
                s_ip_t0 = self.V[t][sname]['s'][self._iassets[ind,t]+1,  self.iexo[ind,t], self._itheta[ind,t]].reshape(nst)
                s_ip_tp = self.V[t][sname]['s'][self._iassets[ind,t]+1,  self.iexo[ind,t], self._itheta[ind,t]+1].reshape(nst)
                
                s_i = (1-wt_p)*s_i_t0 + wt_p*s_i_tp
                s_ip = (1-wt_p)*s_ip_t0 + wt_p*s_ip_tp
                
                
            w_i = self._wnassets[ind,t].reshape(nst)    
            anext_val =  (1-w_i)*s_i + w_i*s_ip 
            
            assert np.all(anext_val >= 0)
            
            self.assets[ind,t+1] = anext_val
            
            # ok what do I do for interpolation now?
            # this is performed when updating the state but do it anyways
            self._iassets[ind,t+1], self._wnassets[ind,t+1] = interp(agrid,anext_val,trim=True)
            
    def iexonext(self,t):
        
        # let's find out new exogenous state
        
        for ist,sname in enumerate(self.state_names):
            is_state = (self.state[:,t]==ist)
            nst = np.sum(is_state)
            
            if nst == 0:
                continue
            
            ind = np.where(is_state)[0]
            sname = self.state_names[ist]
            
            mat = self.setup.exo_mats[sname][t]
            
            iexo_now = self.iexo[ind,t].reshape(nst)
            
            iexo_next = mc_simulate(iexo_now,mat,shocks=None) # import + add shocks     
            self.iexo[ind,t+1] = iexo_next
            
            assert np.all(iexo_next<self.setup.pars['nexo'])
            
    def statenext(self,t):
        
        setup = self.setup
        agrid = setup.agrid
        thetagrid = setup.thetagrid
        
        for ist,sname in enumerate(self.state_names):
            is_state = (self.state[:,t]==ist)
            
            if self.verbose: print('At t = {} count of {} is {}'.format(t,sname,np.sum(is_state)))
            
            if not np.any(is_state):
                continue
            
            
            #Vnext = dict().fromkeys(self.state_names,None)        
            ind = np.where(is_state)[0]
            
            if sname == "Female, single":
                
                # meet a partner
                pmat = self.setup.part_mats['Female, single'][t]
                ic_out = mc_simulate(self.iexo[ind,t],pmat,shocks=None)
                sf = self.assets[ind,t+1] # note that timing is slightly inconsistent
                # we use iexo from t and savings from t+1
                # TODO: fix the seed
                mult_a = np.exp(np.random.normal()*setup.pars['sig_partner_a'])
                sm = mult_a*sf
                iall, izf, izm, ipsi = setup.all_indices(ic_out)
                
                vf, vm, ismar, tht, _ = v_mar(setup,self.V[t+1],sf,sm,iall,combine=False,return_all=True)
                
                
                i_disagree = np.isnan(tht)
                
                i_agree = ~np.array(i_disagree)
                
                print('{} agreed, {} disagreed'.format(np.sum(i_agree),np.sum(i_disagree)))

                
                assert np.all(ismar==i_agree)
                
                if np.any(i_agree):
                    
                    self.theta[ind[i_agree],t+1] = tht[i_agree] # mad indexing
                    self.assets[ind[i_agree],t+1] = sf[i_agree] + sm[i_agree]
                    self.iexo[ind[i_agree],t+1] = iall[i_agree]
                    self.state[ind[i_agree],t+1] = self.state_codes['Couple']
                    
                if np.any(i_disagree):
                    self.theta[ind[i_disagree],t+1] = np.nan
                    # do not touch assets
                    self.iexo[ind[i_disagree],t+1] = izf[i_disagree]
                    self.state[ind[i_disagree],t+1] = self.state_codes['Female, single']
                    
                self._iassets[ind,t+1], self._wnassets[ind,t+1] = interp(agrid,self.assets[ind,t+1],trim=True)
                
                if np.any(i_agree):
                    self._itheta[ind[i_agree],t+1], self._wntheta[ind[i_agree],t+1] = interp(thetagrid,self.theta[ind[i_agree],t+1],trim=True)
                   
                
            elif sname == "Couple":
                
                # by default keep the same theta and weights
                self.theta[ind,t+1] = self.theta[ind,t]
                self._itheta[ind,t+1], self._wntheta[ind,t+1] = self._itheta[ind,t], self._wntheta[ind,t]
                
                # initiate renegotiation
                sc = self.assets[ind,t+1]
                iall, izf, izm, ipsi = self.setup.all_indices(self.iexo[ind,t+1])
                
                v, vf, vm, thetaout, tht_fem, tht_mal = v_ren(setup,self.V[t+1],sc=sc,ind_or_inds=iall,combine=False,return_all=True)
                #print((tht_fem,tht_mal))
                
                # same size as ind
                i_div = np.array(tht_fem > tht_mal)
                i_ren = np.array(~(i_div) & ((self.theta[ind,t] < tht_fem) | (self.theta[ind,t] > tht_mal)))
                i_renf = np.array(i_ren & (self.theta[ind,t] < tht_fem))
                i_renm = np.array(i_ren & (self.theta[ind,t] > tht_mal))
                i_sq  = np.array(~(i_div) & ~(i_ren))
                
                print('{} divorce, {} ren-f, {} ren-m, {} sq'.format(np.sum(i_div),np.sum(i_renf),np.sum(i_renm),np.sum(i_sq)))
                
                
                if np.any(i_div):
                    # TODO: this should replicate the divorce protocol
                    self.assets[ind[i_div],t+1] = 0.5*sc[i_div] 
                    self.theta[ind[i_div],t+1], self._wntheta[ind[i_div],t+1], self._itheta[ind[i_div],t+1] = np.nan, 0.0, 0
                    self.iexo[ind[i_div],t+1] = izf[i_div]
                    self.state[ind[i_div],t+1] = self.state_codes['Female, single']
                
                if np.any(i_ren):
                    self.assets[ind[i_ren],t+1] = sc[i_ren]
                    if np.any(i_renf): self.theta[ind[i_renf],t+1] = tht_fem[i_renf]
                    if np.any(i_renm): self.theta[ind[i_renm],t+1] = tht_fem[i_renm]
                    self.state[ind[i_ren],t+1] = self.state_codes['Couple']
                    
                if np.any(i_sq):
                    self.state[ind[i_sq],t+1] = self.state_codes['Couple']
                    
                self._itheta[ind[i_ren],t+1], self._itheta[ind[i_ren],t+1] = interp(thetagrid,self.theta[ind[i_ren],t+1],trim=True)
                self._iassets[ind[i_div],t+1], self._wnassets[ind[i_div],t+1] = interp(agrid,self.assets[ind[i_div],t+1],trim=True)
                    
                
                #self.state[ind,t+1] = self.state_codes['Couple'] # by default keep couples
                
                
                
                
                
            else:
                raise Exception('unsupported state?')
        
        assert not np.any(np.isnan(self.state[:,t+1]))
                
            
    
        
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
    
    mdl = Model(iterator_name='default-timed',
                divorce_costs={'unilateral_divorce':True})
    mdl.solve_sim()
    mdl.time_statistics()
    