#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This contains things relevant for simulations
"""

import numpy as np


from interp_np import interp            
from mc_tools import mc_simulate
from ren_mar import v_mar, v_ren
from gridvec import VecOnGrid

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
        
        
        
        # initialize assets
        self.gassets = [VecOnGrid(self.setup.agrid,np.zeros(N,dtype=np.float32),trim=True)
                            for _ in range(T)] 
        
        
        # initialize theta
        self.gtheta = [VecOnGrid(self.setup.thetagrid,np.zeros(N,dtype=np.float32),trim=True)
                            for _ in range(T)]
       
        self.theta = np.full((N,T),np.nan,dtype=np.float32)
        self._itheta  = np.zeros((N,T),dtype=np.int32)           
        self._wntheta  = np.zeros((N,T),dtype=np.float32)
        
        
        # initialize iexo
        self.iexo = np.zeros((N,T),np.int32)
        self.iexo[:,0] = np.random.randint(0,self.setup.pars['n_zf'],size=N) # initialize iexo
        
        
        # initialize state
        self.state = np.zeros((N,T),dtype=np.int32)       
        self.state[:,0] = 0  # everyone starts as female
        
        self.state_codes = dict()
        self.has_theta = list()
        for i, name in enumerate(self.setup.state_names):
            self.state_codes[name] = i
            self.has_theta.append((name=='Couple'))
        
        
        
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
            
            
            
            if not use_theta:
                
                anext_val = self.gassets[t].apply( self.V[t][sname]['s'], take = [(1,self.iexo[ind,t])], 
                            pick = ind[0], reshape_i = False).squeeze()
                
            else:
                
                # at this theta
                s_t0 = self.gassets[t].apply( self.V[t][sname]['s'], take = [(1,self.iexo[ind,t]),(2,self._itheta[ind,t])], 
                            pick = ind[0], reshape_i = False).squeeze()
                
                # at next theta
                s_tp = self.gassets[t].apply( self.V[t][sname]['s'], take = [(1,self.iexo[ind,t]),(2,self._itheta[ind,t]+1)], 
                            pick = ind[0], reshape_i = False).squeeze()
                # wegiht of the next theta
                wt_p = self._wntheta[ind,t].reshape(nst)                
                anext_val = (1-wt_p)*s_t0 + wt_p*s_tp                
                
                #tht = VecOnGrid(self.setup.thetagrid, self.theta[:,t],trim=True)                
                
                #q3 = tht.apply(s2,axis=1,pick=ind[0],take=(0, list(range(s2.shape[0])) ),reshape_i=False).squeeze()
                #thetagrid = self.setup.thetagrid
                
                #thetacheck = (thetagrid[self._itheta[ind,t]]*(1-self._wntheta[ind,t]) + \
                #               thetagrid[self._itheta[ind,t]+1]*self._wntheta[ind,t]).squeeze()
                
                #iseq = np.isclose(thetacheck,self.theta[ind,t].squeeze())
                
            assert np.all(anext_val >= 0)
            
            self.gassets[t+1].update(ind[0],anext_val) 
            
            
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
                sf = self.gassets[t+1].val[ind] # note that timing is slightly inconsistent
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
                    self._itheta[ind[i_agree],t+1], self._wntheta[ind[i_agree],t+1] = interp(thetagrid,tht[i_agree],trim=True)
                    
                    self.iexo[ind[i_agree],t+1] = iall[i_agree]
                    self.state[ind[i_agree],t+1] = self.state_codes['Couple']
                    
                    self.gassets[t+1].update(ind[i_agree],sf[i_agree] + sm[i_agree])
                    
                if np.any(i_disagree):
                    #self.theta[ind[i_disagree],t+1] = np.nan
                    #_self._itheta[ind[i_disagree],t+1], _self._wntheta[ind[i_disagree],t+1] = 0, 0.0
                    # do not touch assets
                    self.iexo[ind[i_disagree],t+1] = izf[i_disagree]
                    self.state[ind[i_disagree],t+1] = self.state_codes['Female, single']
                    
                
            elif sname == "Couple":
                
                # by default keep the same theta and weights
                self.theta[ind,t+1] = self.theta[ind,t]
                self._itheta[ind,t+1], self._wntheta[ind,t+1] = self._itheta[ind,t], self._wntheta[ind,t]
                
                # initiate renegotiation
                sc = self.gassets[t+1].val[ind]
                iall, izf, izm, ipsi = self.setup.all_indices(self.iexo[ind,t+1])
                
                v, vf, vm, thetaout, tht_fem, tht_mal = v_ren(setup,self.V[t+1],sc=sc,ind_or_inds=iall,combine=False,return_all=True)
                #print((tht_fem,tht_mal))
                
                # same size as ind
                i_div = np.array(tht_fem > tht_mal)
                i_ren = np.array(~(i_div) & ((self.theta[ind,t] < tht_fem) | (self.theta[ind,t] > tht_mal)))
                i_renf = np.array(i_ren & (self.theta[ind,t] < tht_fem))
                i_renm = np.array(i_ren & (self.theta[ind,t] > tht_mal))
                i_sq  = np.array(~(i_div) & ~(i_ren))
                
                assert np.all(i_ren == ((i_renf) | (i_renm)))
                print('{} divorce, {} ren-f, {} ren-m, {} sq'.format(np.sum(i_div),np.sum(i_renf),np.sum(i_renm),np.sum(i_sq)))
                
                
                if np.any(i_div):
                    # TODO: this should replicate the divorce protocol
                    self.gassets[t+1].update(ind[i_div],0.5*sc[i_div])                    
                    self.theta[ind[i_div],t+1], self._wntheta[ind[i_div],t+1], self._itheta[ind[i_div],t+1] = np.nan, 0.0, 0
                    self.iexo[ind[i_div],t+1] = izf[i_div]
                    self.state[ind[i_div],t+1] = self.state_codes['Female, single']
                
                
                if np.any(i_ren):
                    if np.any(i_renf): 
                        self.theta[ind[i_renf],t+1] = tht_fem[i_renf]
                        self._itheta[ind[i_renf],t+1], self._wntheta[ind[i_renf],t+1] = interp(thetagrid,self.theta[ind[i_renf],t+1],trim=True)
                    if np.any(i_renm):
                        self.theta[ind[i_renm],t+1] = tht_mal[i_renm]
                        self._itheta[ind[i_renm],t+1], self._wntheta[ind[i_renm],t+1] = interp(thetagrid,self.theta[ind[i_renm],t+1],trim=True)
                    self.state[ind[i_ren],t+1] = self.state_codes['Couple']
                    
                
                if np.any(i_sq):
                    self.state[ind[i_sq],t+1] = self.state_codes['Couple']
                    
                
            
            else:
                raise Exception('unsupported state?')
        
        assert not np.any(np.isnan(self.state[:,t+1]))
                