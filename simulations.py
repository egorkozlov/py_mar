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
        
        
        self.gassets = [VecOnGrid(self.setup.agrid,0.01*np.ones(N,dtype=np.float32),trim=True)
                            for _ in range(T)]
        
        self._iassets = np.zeros((N,T),np.int32)            
        self._wnassets = np.zeros((N,T),np.float32)
        
        self._itheta  = np.zeros((N,T),dtype=np.int32)           
        self._wntheta  = np.zeros((N,T),dtype=np.float32)
        
        
        self.assets = 0.01*np.ones((N,T),dtype=np.float32)
        self._iassets[:,0], self._wnassets[:,0] = \
            interp(self.setup.agrid,self.assets[:,0],return_wnext=True,trim=True)
        
        
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
            
            
            w_i = self._wnassets[ind,t].reshape(nst)  
            
            if not use_theta:
                s_i  = self.V[t][sname]['s'][self._iassets[ind,t],  self.iexo[ind,t]].reshape(nst)
                s_ip = self.V[t][sname]['s'][self._iassets[ind,t]+1,self.iexo[ind,t]].reshape(nst)
                
                anext_val =  (1-w_i)*s_i + w_i*s_ip 
                
                
                
                q = self.gassets[t].apply( self.V[t][sname]['s'], take = [(1,self.iexo[ind,t])], 
                            pick = ind[0], reshape_i = False).squeeze()
                
                assert np.all( np.abs(anext_val-q) < 1e-5 )
                
                
                
            else:
                wt_p = self._wntheta[ind,t].reshape(nst)                
                # bilinear interpolation                
                s_i_t0 = self.V[t][sname]['s'][self._iassets[ind,t],  self.iexo[ind,t], self._itheta[ind,t]].reshape(nst)
                s_i_tp = self.V[t][sname]['s'][self._iassets[ind,t],  self.iexo[ind,t], self._itheta[ind,t]+1].reshape(nst)
                s_ip_t0 = self.V[t][sname]['s'][self._iassets[ind,t]+1,  self.iexo[ind,t], self._itheta[ind,t]].reshape(nst)
                s_ip_tp = self.V[t][sname]['s'][self._iassets[ind,t]+1,  self.iexo[ind,t], self._itheta[ind,t]+1].reshape(nst)
                
                
                s_i = (1-wt_p)*s_i_t0 + wt_p*s_i_tp
                s_ip = (1-wt_p)*s_ip_t0 + wt_p*s_ip_tp
                
                anext_val =  (1-w_i)*s_i + w_i*s_ip 
                
                
                aq0 = agrid[self._iassets[ind,t]]*(1-w_i) + agrid[self._iassets[ind,t]+1]*w_i
                aq1 = self.gassets[t].val[ind]
                
                try:
                    assert np.all(np.abs( aq0 -  aq1  ) < 1e-5)
                except:
                    print((t,sname))
                    #print(aq0)
                    print((aq0-aq1))
                    #print(aq1)
                    
                    assert False, 'aq are not equal'
                
                s_t0 = self.gassets[t].apply( self.V[t][sname]['s'], take = [(1,self.iexo[ind,t]),(2,self._itheta[ind,t])], 
                            pick = ind[0], reshape_i = False).squeeze()
                
                assert np.all( np.abs(  s_t0 - (1-w_i)*s_i_t0 - w_i*s_ip_t0  ) < 1e-5 )
                
                s_tp = self.gassets[t].apply( self.V[t][sname]['s'], take = [(1,self.iexo[ind,t]),(2,self._itheta[ind,t]+1)], 
                            pick = ind[0], reshape_i = False).squeeze()
                
                assert np.all( np.abs(  s_tp - (1-w_i)*s_i_tp - w_i*s_ip_tp  ) < 1e-5 )
                
                
                s2 = self.gassets[t].apply( self.V[t][sname]['s'], take = [(1,self.iexo[ind,t])], 
                            pick = ind[0], reshape_i = False).squeeze()
                
                # try to do different indexing schemes
                s2_t0 = s2[list(range(s2.shape[0])), self._itheta[ind,t].squeeze()]                
                assert np.all( np.abs( s_t0 - s2_t0) < 1e-5 )                
                s2_tp = s2[list(range(s2.shape[0])),self._itheta[ind,t].squeeze()+1]                
                assert np.all( np.abs( s_tp - s2_tp) < 1e-5 )
                
                
                
                
                
                q = (1-wt_p)*s_t0 + wt_p*s_tp                
                tht = VecOnGrid(self.setup.thetagrid, self.theta[:,t],trim=True)                
                
                q3 = tht.apply(s2,axis=1,pick=ind[0],take=(0, list(range(s2.shape[0])) ),reshape_i=False).squeeze()
                thetagrid = self.setup.thetagrid
                
                thetacheck = (thetagrid[self._itheta[ind,t]]*(1-self._wntheta[ind,t]) + \
                               thetagrid[self._itheta[ind,t]+1]*self._wntheta[ind,t]).squeeze()
                iseq = ( np.abs(thetacheck - self.theta[ind,t].squeeze()) < 1e-5 )
                                    
                                 
                try:
                    assert np.all(iseq)
                except:
                    print((t,sname))
                    print(thetacheck.shape)
                    print(thetacheck[np.where(~(iseq))])
                    print(self.theta[ind,t].squeeze()[np.where(~(iseq))])
                    assert False
                    
                    
                
                q2 = ( (1-wt_p)*s2_t0 + wt_p*s2_tp ).squeeze()
                
                try:
                    assert np.all( np.abs(q3-q2) < 1e-5 )
                except:
                    print('q3 is not q2')
                    print(q3-q2)
                    print((q3.shape,q2.shape))                    
                    assert False
                    
                    
                try:
                    assert np.all( np.abs(anext_val - q2) < 1e-5)
                except:
                    #print(anext_val)
                    print('anext is not q2')
                    print(anext_val - q2)
                    assert False
                    
                
                
                
                
              
            
            assert np.all(anext_val >= 0)
            
            self.assets[ind,t+1] = anext_val
            self.gassets[t+1].update(ind[0],anext_val) 
            
            
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
                    self._itheta[ind[i_agree],t+1], self._wntheta[ind[i_agree],t+1] = interp(thetagrid,tht[i_agree],trim=True)
                    
                    self.assets[ind[i_agree],t+1] = sf[i_agree] + sm[i_agree]
                    self._iassets[ind[i_agree],t+1], self._wnassets[ind[i_agree],t+1] = interp(agrid,self.assets[ind[i_agree],t+1],trim=True)
                
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
                
                
                assert np.all(i_ren == ((i_renf) | (i_renm)))
                
                print('{} divorce, {} ren-f, {} ren-m, {} sq'.format(np.sum(i_div),np.sum(i_renf),np.sum(i_renm),np.sum(i_sq)))
                
                
                if np.any(i_div):
                    # TODO: this should replicate the divorce protocol
                    self.assets[ind[i_div],t+1] = 0.5*sc[i_div] 
                    self.theta[ind[i_div],t+1], self._wntheta[ind[i_div],t+1], self._itheta[ind[i_div],t+1] = np.nan, 0.0, 0
                    self.iexo[ind[i_div],t+1] = izf[i_div]
                    self.state[ind[i_div],t+1] = self.state_codes['Female, single']
                
                    self.gassets[t+1].update(ind[i_div],0.5*sc[i_div])
                    
                if np.any(i_ren):
                    
                    self.assets[ind[i_ren],t+1] = sc[i_ren]
                    if np.any(i_renf): 
                        self.theta[ind[i_renf],t+1] = tht_fem[i_renf]
                        self._itheta[ind[i_renf],t+1], self._wntheta[ind[i_renf],t+1] = interp(thetagrid,self.theta[ind[i_renf],t+1],trim=True)
                    if np.any(i_renm):
                        self.theta[ind[i_renm],t+1] = tht_mal[i_renm]
                        self._itheta[ind[i_renm],t+1], self._wntheta[ind[i_renm],t+1] = interp(thetagrid,self.theta[ind[i_renm],t+1],trim=True)
                    
                    
                    self.state[ind[i_ren],t+1] = self.state_codes['Couple']
                    
                if np.any(i_sq):
                    self.state[ind[i_sq],t+1] = self.state_codes['Couple']
                    
                
                self._iassets[ind[i_div],t+1], self._wnassets[ind[i_div],t+1] = interp(agrid,self.assets[ind[i_div],t+1],trim=True)
                    
                
                
                
            else:
                raise Exception('unsupported state?')
        
        assert not np.any(np.isnan(self.state[:,t+1]))
                