#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This contains things relevant for simulations
"""

import numpy as np


from interp_np import interp            
from mc_tools import mc_simulate
from ren_mar import v_mar2, v_ren2
from gridvec import VecOnGrid

class Agents:
    
    
    def __init__(self,M,N=1000,T=None,verbose=True):
        if T is None:
            T = M.setup.pars['T']
            
            
        np.random.seed(18) # TODO: this should be replaced by explicitly supplying shocks  
            
        # take the stuff from the model and arguments
        # note that this does not induce any copying just creates links
        self.M = M
        self.V = M.V
        self.setup = M.setup
        self.state_names = self.setup.state_names
        self.N = N
        self.T = T
        self.verbose = verbose
        self.timer = M.time
        
        
        
        # initialize assets
        self.gassets = [VecOnGrid(self.setup.agrid,np.zeros(N,dtype=np.float32),trim=True)
                            for _ in range(T)] 
        # initialize theta
        self.gtheta = [VecOnGrid(self.setup.thetagrid,-1*np.ones(N,dtype=np.float32),trim=True)
                            for _ in range(T)]
        
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
            self.has_theta.append((name=='Couple, C' or name=='Couple, M'))
        
        
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
            
            ind = np.where(is_state)[0]
            
            if not use_theta:
                
                anext = self.gassets[t].apply( self.V[t][sname]['s'], take = [(1,self.iexo[ind,t])], 
                            pick = ind, reshape_i = False)
                
            else:
                
                # interpolate in both assets and theta
                # function apply_2dim is experimental but I checked it at this setup
                
                anext = self.gassets[t].apply_2dim(self.V[t][sname]['s'],
                                                     apply_first=self.gtheta[t],
                                                     axis_first=2,
                                                     axis_this=0,
                                                     take = [(1,self.iexo[ind,t])],
                                                     pick = ind,
                                                     reshape_i = False
                                                    )
                
                
                '''
                # This is a consistency check about what if we interpolate
                # in different order
                
                anext_alt = self.gtheta[t].apply_2dim(self.V[t][sname]['s'],
                                                     apply_first=self.gassets[t],
                                                     axis_first=0,
                                                     axis_this=2,
                                                     take = [(1,self.iexo[ind,t])],
                                                     pick = ind,
                                                     reshape_i = False
                                                    )
                
                
                assert np.allclose(anext,anext_alt)
                
                '''
                
                
                
            assert np.all(anext >= 0)
            
            self.gassets[t+1].update(ind,anext) 
            
            
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
        
        for ist,sname in enumerate(self.state_names):
            is_state = (self.state[:,t]==ist)
            
            if self.verbose: print('At t = {} count of {} is {}'.format(t,sname,np.sum(is_state)))
            
            if not np.any(is_state):
                continue
            
            ind = np.where(is_state)[0]
            
            nind = ind.size
            
            
            
            if sname == "Female, single":
                # TODO: this is temporary version, it computes partners for
                # everyone and after that imposes meet / no meet, this should
                # not happen.
                
                # meet a partner
                pmat = self.setup.part_mats['Female, single'][t]
                pmeet = self.setup.pars['pmeet']
                ic_out = mc_simulate(self.iexo[ind,t],pmat,shocks=None)
                iall, izf, izm, ipsi = setup.all_indices(ic_out)
                
                sf = self.gassets[t+1].val[ind] # note that timing is slightly inconsistent
                # we use iexo from t and savings from t+1
                # TODO: fix the seed
                
                mult_a = np.exp(np.random.normal()*setup.pars['sig_partner_a'])
                sm = mult_a*sf
                
                # compute for everyone
                
                
                # cohabitation
                (vf_c, vm_c, iscoh, tht_c, _), nbs_c = v_mar2(setup,self.V[t+1],False,sf,sm,iall,combine=False,return_all=True)
                # marriage
                (vf_m, vm_m, ismar, tht_m, _), nbs_m = v_mar2(setup,self.V[t+1],True,sf,sm,iall,combine=False,return_all=True)
                
                
                i_nomeet =  np.array( np.random.rand(nind) > pmeet )
                
                
                i_pot_mar = ~np.isnan(tht_m)
                
                assert np.all( i_pot_mar == (nbs_m>0) )
                
                i_pot_coh = ~np.isnan(tht_c)
                
                assert np.all( i_pot_coh == (nbs_c > 0) )
                
                
                i_disagree = (~i_pot_mar) & (~i_pot_coh)
                i_disagree_or_nomeet = (i_disagree) | (i_nomeet)
                i_agree = ~np.array(i_disagree_or_nomeet)
                
                i_agree_mar = (i_agree) & (nbs_m>=nbs_c) 
                i_agree_coh = (i_agree) & (nbs_c> nbs_m)
                
                assert np.all(~i_nomeet[i_agree])
                
                
                nmar, ncoh, ndis, nnom = np.sum(i_agree_mar),np.sum(i_agree_coh),np.sum(i_disagree_or_nomeet),np.sum(i_nomeet)
                ntot = sum((nmar, ncoh, ndis, nnom))
                
                print('{} mar, {} coh,  {} disagreed, {} did not meet ({} total)'.format(nmar,ncoh,ndis,nnom,ntot))
                #assert np.all(ismar==(i_agree )
                
                if np.any(i_agree_mar):
                    
                    self.gtheta[t+1].update(ind[i_agree_mar],tht_m[i_agree_mar])
                    self.iexo[ind[i_agree_mar],t+1] = iall[i_agree_mar]
                    self.state[ind[i_agree_mar],t+1] = self.state_codes['Couple, M']
                    self.gassets[t+1].update(ind[i_agree_mar],sf[i_agree_mar] + sm[i_agree_mar])
                    
                if np.any(i_agree_coh):
                    
                    self.gtheta[t+1].update(ind[i_agree_coh],tht_c[i_agree_coh])
                    self.iexo[ind[i_agree_coh],t+1] = iall[i_agree_coh]
                    self.state[ind[i_agree_coh],t+1] = self.state_codes['Couple, C']
                    self.gassets[t+1].update(ind[i_agree_coh],sf[i_agree_coh] + sm[i_agree_coh])
                    
                
                    
                if np.any(i_disagree_or_nomeet):
                    # do not touch assets
                    self.iexo[ind[i_disagree_or_nomeet],t+1] = izf[i_disagree_or_nomeet]
                    self.state[ind[i_disagree_or_nomeet],t+1] = self.state_codes['Female, single']
                    
            elif sname == "Couple, M" or sname == "Couple, C":
                
                # by default keep the same theta and weights
                thetanow = self.gtheta[t].val[ind]
                self.gtheta[t+1].update(ind,thetanow)
                
                # initiate renegotiation
                sc = self.gassets[t+1].val[ind]
                iall, izf, izm, ipsi = self.setup.all_indices(self.iexo[ind,t+1])
                
                v, vf, vm, thetaout, tht_fem, tht_mal = v_ren2(setup,self.V[t+1],True,t,sc=sc,ind_or_inds=iall,combine=False,return_all=True)
                
                
                # same size as ind
                # categorizes couples
                i_div = np.array(tht_fem > tht_mal)
                i_ren = np.array(~(i_div) & ((thetanow < tht_fem) | (thetanow > tht_mal)))
                i_renf = np.array(i_ren & (thetanow < tht_fem))
                i_renm = np.array(i_ren & (thetanow > tht_mal))
                i_sq  = np.array(~(i_div) & ~(i_ren))
                
                assert np.all(i_ren == ((i_renf) | (i_renm)))
                print('{} divorce, {} ren-f, {} ren-m, {} sq'.format(np.sum(i_div),np.sum(i_renf),np.sum(i_renm),np.sum(i_sq))                     )
                
                
                
                zf_grid = self.setup.exo_grids['Female, single'][t]
                zm_grid = self.setup.exo_grids['Male, single'][t]
                
                
                
                if np.any(i_div):
                    
                    income_fem = np.exp(zf_grid[izf[i_div]])
                    income_mal = np.exp(zm_grid[izm[i_div]])
                    
                    income_share_fem = income_fem / (income_fem + income_mal)
                    
                    costs = self.setup.div_costs if sname == 'Couple, M' else self.setup.sep_costs
                               
                    share_f, share_m = costs.shares_if_split(income_share_fem)
                    
                    sf = share_f*sc[i_div]
                    
                    self.gassets[t+1].update(ind[i_div],sf)                    
                    self.gtheta[t+1].update(ind[i_div],-1.0)
                    self.iexo[ind[i_div],t+1] = izf[i_div]
                    self.state[ind[i_div],t+1] = self.state_codes['Female, single']
                
                if np.any(i_ren):
                    if np.any(i_renf): 
                        self.gtheta[t+1].update(ind[i_renf],tht_fem[i_renf])
                    if np.any(i_renm):
                        self.gtheta[t+1].update(ind[i_renm],tht_mal[i_renm])
                    self.state[ind[i_ren],t+1] = self.state_codes[sname]
                    
                if np.any(i_sq):
                    self.state[ind[i_sq],t+1] = self.state_codes[sname]
                    # do not touch theta as already updated
                
            
            else:
                raise Exception('unsupported state?')
        
        assert not np.any(np.isnan(self.state[:,t+1]))
                