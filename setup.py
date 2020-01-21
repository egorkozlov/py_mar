#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This contains things relevant for setting up the model
"""

import numpy as np

from rw_approximations import rouw_nonst
from mc_tools import combine_matrices_two_lists, int_prob,cut_matrix
from scipy.stats import norm
from collections import namedtuple
from gridvec import VecOnGrid

from scipy import sparse



class ModelSetup(object):
    def __init__(self,nogrid=False,divorce_costs='Default',separation_costs='Default',**kwargs): 
        p = dict()       
        period_year=2#this can be 1,2,3 or 6
        T = int(52/period_year)
        Tret = int(42/period_year) # first period when the agent is retired
        Tbef=int(2/period_year)
        p['py']=period_year
        p['T'] = T
        p['Tret'] = Tret
        p['Tbef'] = Tbef
        p['sig_zf_0']  = 0.25
        p['sig_zf']    = 0.075
        p['n_zf_t']      = [7]*Tret + [1]*(T-Tret)
        p['sig_zm_0']  = 0.25
        p['sig_zm']    = 0.05
        p['n_zm_t']      = [5]*Tret + [1]*(T-Tret)
        p['sigma_psi_init'] = 0.28
        p['sigma_psi']   = 0.11
        p['R_t'] = [1.04**period_year]*T
        p['n_psi_t']     = [12]*T
        p['beta_t'] = [0.98**period_year]*T
        p['A'] = 1.0 # consumption in couple: c = (1/A)*[c_f^(1+rho) + c_m^(1+rho)]^(1/(1+rho))
        p['crra_power'] = 2
        p['couple_rts'] = 0.00
        p['sig_partner_a'] = 0.1
        p['sig_partner_z'] = 0.4
        p['m_bargaining_weight'] = 0.5
        p['pmeet'] = 0.5
        
        p['target_x_share'] = 0.6
        p['target_x_theta'] = 0.4
        
        
        
        p['wret'] = 0.8
        p['uls'] = 0.3
        #p['uls_x'] = 0.0
        p['pls'] = 0.8
        
        
        p['u_shift_mar'] = 0.0
        p['u_shift_coh'] = 0.0
        
        
        p['f_wage_trend'] = [0.0 + 0.0*min(t,Tret) - 0.000*(min(t,Tret)**2) for t in range(T)]
        p['m_wage_trend'] = [0.0 + 0.0*min(t,Tret) - 0.000*(min(t,Tret)**2) for t in range(T)]
        
        
        
        for key, value in kwargs.items():
            assert (key in p), 'wrong name?'
            p[key] = value
        
        #Get the probability of meeting, adjusting for year-period
        p_meet=p['pmeet']
        for j in range(period_year-1):
            p_meet=p_meet+(1-p_meet)*p['pmeet']
            
        p['pmeet_t'] = [p_meet]*T
        
        
        self.pars = p # !!! important
        
        
        
        shr = p['target_x_share']
        tht = p['target_x_theta']
        A_theta = self.u_mult_just_c(np.array([tht]))
        
        ulsx = float( A_theta * (shr/(1-shr))**(p['crra_power']) )
        self.pars['uls_x'] = ulsx
        p['uls_x'] = ulsx
        
        
        self.dtype = np.float32 # type for all floats
        
       
        
        # relevant for integration
        self.state_names = ['Female, single','Male, single','Couple, M', 'Couple, C']
        
        # female labor supply
        self.ls_levels = np.array([0.2,1.0],dtype=self.dtype)
        self.ls_u0 = np.array([p['uls'],0.0],dtype=self.dtype)
        self.ls_pdown = np.array([p['pls'],0.0],dtype=self.dtype)
        self.nls = len(self.ls_levels)     
        
        
        
        
        
        self.ls_u1 = np.array([p['uls_x'],p['uls_x']],dtype=self.dtype)
        
        
        
        
        
        #Cost of Divorce
        if divorce_costs == 'Default':
            # by default the costs are set in the bottom
            self.div_costs = DivorceCosts()
        else:
            if isinstance(divorce_costs,dict):
                # you can feed in arguments to DivorceCosts
                self.div_costs = DivorceCosts(**divorce_costs)
            else:
                # or just the output of DivorceCosts
                assert isinstance(divorce_costs,DivorceCosts)
                self.div_costs = divorce_costs
                
        #Cost of Separation
        if separation_costs == 'Default':
            # by default the costs are set in the bottom
            self.sep_costs = DivorceCosts()
        else:
            if isinstance(separation_costs,dict):
                # you can feed in arguments to DivorceCosts
                self.sep_costs = DivorceCosts(**divorce_costs)
            else:
                # or just the output of DivorceCosts
                assert isinstance(separation_costs,DivorceCosts)
                self.sep_costs = separation_costs
            
        # exogrid should be deprecated
        if not nogrid:
        
            exogrid = dict()
            
            
            # let's approximate three Markov chains
            # this sets up exogenous grid
            
            # FIXME: this uses number of points from 0th entry. 
            # in principle we can generalize this
            
            exogrid['zf_t'],  exogrid['zf_t_mat'] = rouw_nonst(p['T'],p['sig_zf']*period_year,p['sig_zf_0'],p['n_zf_t'][0])
            exogrid['zm_t'],  exogrid['zm_t_mat'] = rouw_nonst(p['T'],p['sig_zm']*period_year,p['sig_zm_0'],p['n_zm_t'][0])
            
            for t in range(Tret,T):
                exogrid['zf_t'][t] = np.array([np.log(p['wret'])])
                exogrid['zm_t'][t] = np.array([np.log(p['wret'])])
                exogrid['zf_t_mat'][t] = np.atleast_2d(1.0)
                exogrid['zm_t_mat'][t] = np.atleast_2d(1.0)
                
                
            # fix transition from non-retired to retired    
            exogrid['zf_t_mat'][Tret-1] = np.ones((p['n_zf_t'][Tret-1],1))
            exogrid['zm_t_mat'][Tret-1] = np.ones((p['n_zm_t'][Tret-1],1))
            
            exogrid['psi_t'], exogrid['psi_t_mat'] = rouw_nonst(p['T'],p['sigma_psi']*period_year,p['sigma_psi_init'],p['n_psi_t'][0])
            
            zfzm, zfzmmat = combine_matrices_two_lists(exogrid['zf_t'], exogrid['zm_t'], exogrid['zf_t_mat'], exogrid['zm_t_mat'])
            all_t, all_t_mat = combine_matrices_two_lists(zfzm,exogrid['psi_t'],zfzmmat,exogrid['psi_t_mat'])
            all_t_mat_sparse_T = [sparse.csc_matrix(D.T) if D is not None else None for D in all_t_mat]
            
            
            #Create a new bad version of transition matrix p(zf_t)
            
            
            zf_bad = [cut_matrix(exogrid['zf_t_mat'][t]) if t < Tret -1 
                          else (exogrid['zf_t_mat'][t] if t < T - 1 else None) 
                              for t in range(self.pars['T'])]
            
            zf_t_mat_down = zf_bad
            
            zfzm, zfzmmat = combine_matrices_two_lists(exogrid['zf_t'], exogrid['zm_t'], zf_t_mat_down, exogrid['zm_t_mat'])
            all_t_down, all_t_mat_down = combine_matrices_two_lists(zfzm,exogrid['psi_t'],zfzmmat,exogrid['psi_t_mat'])
            all_t_mat_down_sparse_T = [sparse.csc_matrix(D.T) if D is not None else None for D in all_t_mat_down]
            
            
            
            all_t_mat_by_l = [ [(1-p)*m + p*md if m is not None else None 
                                for m , md in zip(all_t_mat,all_t_mat_down)]
                               for p in self.ls_pdown ]
            
            all_t_mat_by_l_spt = [ [(1-p)*m + p*md if m is not None else None
                                    for m, md in zip(all_t_mat_sparse_T,all_t_mat_down_sparse_T)]
                               for p in self.ls_pdown ]
            
            
            
            exogrid['all_t_mat_by_l'] = all_t_mat_by_l
            exogrid['all_t_mat_by_l_spt'] = all_t_mat_by_l_spt
            
            exogrid['all_t'] = all_t
            
            Exogrid_nt = namedtuple('Exogrid_nt',exogrid.keys())
            
            self.exogrid = Exogrid_nt(**exogrid)
            self.pars['nexo_t'] = [v.shape[0] for v in all_t]
            
            #assert False
            
        #Grid Couple
        self.na = 40
        self.amin = 0
        self.amax = 80
        #self.agrid_c = np.linspace(self.amin,self.amax,self.na,dtype=self.dtype)
        tune=1.5
        self.agrid_c = np.geomspace(self.amin+tune,self.amax+tune,num=self.na)-tune
        
        # this builds finer grid for potential savings
        s_between = 10 # default numer of points between poitns on agrid
        s_da_min = 0.001 # minimal step (does not create more points)
        s_da_max = 0.1 # maximal step (creates more if not enough)
        
        self.sgrid_c = build_s_grid(self.agrid_c,s_between,s_da_min,s_da_max)
        self.vsgrid_c = VecOnGrid(self.agrid_c,self.sgrid_c)
        
        #Grid Single
        self.amin_s = 0
        self.amax_s = self.amax/2.0
        #self.agrid_s = np.linspace(self.amin_s,self.amax_s,self.na,dtype=self.dtype)
        tune_s=1.5
        self.agrid_s = np.geomspace(self.amin_s+tune_s,self.amax_s+tune_s,num=self.na)-tune_s
        
        self.sgrid_s = build_s_grid(self.agrid_s,s_between,s_da_min,s_da_max)
        self.vsgrid_s = VecOnGrid(self.agrid_s,self.sgrid_s)
        
        # grid for theta
        self.ntheta = 11
        self.thetamin = 0.01
        self.thetamax = 0.99
        self.thetagrid = np.linspace(self.thetamin,self.thetamax,self.ntheta,dtype=self.dtype)
        
        
        
        
        
        
        # construct finer grid for bargaining
        ntheta_fine = 10*self.ntheta # actual number may be a bit bigger
        self.thetagrid_fine = np.unique(np.concatenate( (self.thetagrid,np.linspace(self.thetamin,self.thetamax,ntheta_fine,dtype=self.dtype)) ))
        self.ntheta_fine = self.thetagrid_fine.size
        
        i_orig = list()
        
        for theta in self.thetagrid:
            assert theta in self.thetagrid_fine
            i_orig.append(np.where(self.thetagrid_fine==theta)[0])
            
        assert len(i_orig) == self.thetagrid.size
        # allows to recover original gird points on the fine grid        
        self.theta_orig_on_fine = np.array(i_orig).flatten()
        self.v_thetagrid_fine = VecOnGrid(self.thetagrid,self.thetagrid_fine)
        # precomputed object for interpolation

            
        
        
        


        self.exo_grids = {'Female, single':exogrid['zf_t'],
                          'Male, single':exogrid['zm_t'],
                          'Couple, M':exogrid['all_t'],
                          'Couple, C':exogrid['all_t']}
        self.exo_mats = {'Female, single':exogrid['zf_t_mat'],
                          'Male, single':exogrid['zm_t_mat'],
                          'Couple, M':exogrid['all_t_mat_by_l'],
                          'Couple, C':exogrid['all_t_mat_by_l']} # sparse version?
        
        
        self.utility_shifters = {'Female, single':0.0,
                                 'Male, single':0.0,
                                 'Couple, M':p['u_shift_mar'],
                                 'Couple, C':p['u_shift_coh']}
        
        
        # this pre-computes transition matrices for meeting a partner
        zf_t_partmat = [self.mar_mats_iexo(t,female=True) if t < p['T'] - 1 else None 
                            for t in range(p['T'])]
        zm_t_partmat = [self.mar_mats_iexo(t,female=False) if t < p['T'] - 1 else None 
                            for t in range(p['T'])]
        
        self.part_mats = {'Female, single':zf_t_partmat,
                          'Male, single':  zm_t_partmat,
                          'Couple, M': None,
                          'Couple, C': None} # last is added for consistency
        
        self.mar_mats_assets()
        
        self.mar_mats_combine()
        
    def mar_mats_assets(self,npoints=4,abar=0.1):
        # for each grid point on single's grid it returns npoints positions
        # on (potential) couple's grid's and assets of potential partner 
        # (that can be off grid) and correpsonding probabilities. 
        
        
        na = self.agrid_s.size
        
        agrid_s = self.agrid_s
        agrid_c = self.agrid_c
        
        s_a_partner = self.pars['sig_partner_a']
        
        
        prob_a_mat = np.zeros((na,npoints),dtype=self.dtype)
        i_a_mat = np.zeros((na,npoints),dtype=np.int16)
        
        
        
        for ia, a in enumerate(agrid_s):
            lagrid_t = np.zeros_like(agrid_c)
            
            i_neg = (agrid_c <= max(abar,a) - 1e-3)
            
            lagrid_t[~i_neg] = np.log(2e-3 + (agrid_c[~i_neg] - a)/max(abar,a))
            lmin = lagrid_t[~i_neg].min()
            # just fill with very negative values so this is never chosen
            lagrid_t[i_neg] = lmin - s_a_partner*10 - \
                s_a_partner*np.flip(np.arange(i_neg.sum())) 
            
            # TODO: this needs to be checked
            p_a = int_prob(lagrid_t,mu=0,sig=s_a_partner,n_points=npoints)
            i_pa = (-p_a).argsort()[:npoints] # this is more robust then nonzero
            p_pa = p_a[i_pa]
            prob_a_mat[ia,:] = p_pa
            i_a_mat[ia,:] = i_pa
        
        
        self.prob_a_mat = prob_a_mat
        self.i_a_mat = i_a_mat
            
            
        
    
    def mar_mats_iexo(self,t,female=True,trim_lvl=0.001):
        # TODO: check timing
        # this returns transition matrix for single agents into possible couples
        # rows are single's states
        # columnts are couple's states
        # you have to transpose it if you want to use it for integration
        setup = self
        
        nexo = setup.pars['nexo_t'][t]
        sigma_psi_init = setup.pars['sigma_psi_init']
        sig_z_partner = setup.pars['sig_partner_z']
        psi_couple = setup.exogrid.psi_t[t+1]
        
        
        if female:
            nz_single = setup.exogrid.zf_t[t].shape[0]
            p_mat = np.empty((nexo,nz_single))
            z_own = setup.exogrid.zf_t[t]
            n_zown = z_own.shape[0]
            z_partner = setup.exogrid.zm_t[t+1]
            zmat_own = setup.exogrid.zf_t_mat[t]
        else:
            nz_single = setup.exogrid.zm_t[t].shape[0]
            p_mat = np.empty((nexo,nz_single))
            z_own = setup.exogrid.zm_t[t]
            n_zown = z_own.shape[0]
            z_partner = setup.exogrid.zf_t[t+1]
            zmat_own = setup.exogrid.zm_t_mat[t]    
            
        def ind_conv(a,b,c): return setup.all_indices(t,(a,b,c))[0]
        
        
        for iz in range(n_zown):
            p_psi = int_prob(psi_couple,mu=0,sig=sigma_psi_init)
            if female:
                p_zm  = int_prob(z_partner, mu=z_own[iz],sig=sig_z_partner)
                p_zf  = zmat_own[iz,:]
            else:
                p_zf  = int_prob(z_partner, mu=z_own[iz],sig=sig_z_partner)
                p_zm  = zmat_own[iz,:]
            #sm = sf
        
            p_vec = np.zeros(nexo)
            
            for izf, p_zf_i in enumerate(p_zf):
                if p_zf_i < trim_lvl: continue
            
                for izm, p_zm_i in enumerate(p_zm):
                    if p_zf_i*p_zm_i < trim_lvl: continue
                
                    for ipsi, p_psi_i in enumerate(p_psi):                    
                        p = p_zf_i*p_zm_i*p_psi_i
                        
                        if p > trim_lvl:
                            p_vec[ind_conv(izf,izm,ipsi)] = p    
                            
            assert np.any(p_vec>trim_lvl), 'Everything is zero?'              
            p_vec = p_vec / np.sum(p_vec)
            p_mat[:,iz] = p_vec
            
        return p_mat.T
    
    
    def mar_mats_combine(self):
        # for time moment and each position in grid for singles (ia,iz)
        # it computes probability distribution over potential matches
        # this is relevant for testing and simulations
        
        
        pmat_a = self.prob_a_mat
        imat_a = self.i_a_mat
        
        self.matches = dict()
        
        for female in [True,False]:
            desc = 'Female, single' if female else 'Male, single'
            
            pmats = self.part_mats[desc] 
            
            
            match_matrix = list()
            
            
            for t in range(self.pars['T']-1):
                pmat_iexo = pmats[t] # nz X nexo
                # note that here we do not use transpose
                
                nz = pmat_iexo.shape[0]
                
                inds = np.where( np.any(pmat_iexo>0,axis=0) )[0]
                
                npos_iexo = inds.size
                npos_a = pmat_a.shape[1]
                npos = npos_iexo*npos_a
                pmatch = np.zeros((self.na,nz,npos),dtype=self.dtype)
                iamatch = np.zeros((self.na,nz,npos),dtype=np.int32)
                iexomatch = np.zeros((self.na,nz,npos),dtype=np.int32)
                
                i_conv = np.zeros((npos_iexo,npos_a),dtype=np.int32)
                
                for ia in range(npos_a):
                    i_conv[:,ia] = np.arange(npos_iexo*ia,npos_iexo*(ia+1))
                 
                
                for iz in range(nz):
                    probs = pmat_iexo[iz,inds]
                    
                    for ia in range(npos_a):
                        
                        pmatch[:,iz,(npos_iexo*ia):(npos_iexo*(ia+1))] = \
                            (pmat_a[:,ia][:,None])*(probs[None,:])
                        iamatch[:,iz,(npos_iexo*ia):(npos_iexo*(ia+1))] = \
                            imat_a[:,ia][:,None]
                        iexomatch[:,iz,(npos_iexo*ia):(npos_iexo*(ia+1))] = \
                            inds[None,:]
                            
                        
                assert np.allclose(np.sum(pmatch,axis=2),1.0)
                match_matrix.append({'p':pmatch,'ia':iamatch,'iexo':iexomatch,'iconv':i_conv})
                    
            self.matches[desc] = match_matrix
         
        
    
    
    def all_indices(self,t,ind_or_inds=None):
        
        # just return ALL indices if no argument is called
        if ind_or_inds is None: 
            ind_or_inds = np.array(range(self.pars['nexo_t'][t]))
        
        if isinstance(ind_or_inds,tuple):
            izf,izm,ipsi = ind_or_inds
            ind = izf*self.pars['n_zm_t'][t]*self.pars['n_psi_t'][t] + izm*self.pars['n_psi_t'][t] + ipsi
        else:
            ind = ind_or_inds
            izf = ind // (self.pars['n_zm_t'][t]*self.pars['n_psi_t'][t])
            izm = (ind - izf*self.pars['n_zm_t'][t]*self.pars['n_psi_t'][t]) // self.pars['n_psi_t'][t]
            ipsi = ind - izf*self.pars['n_zm_t'][t]*self.pars['n_psi_t'][t] - izm*self.pars['n_psi_t'][t]
            
        return ind, izf, izm, ipsi

    
    # functions u_mult and c_mult are meant to be shape-perservings
    
    # this is multiplier for utility of just-consumption
    def u_mult_just_c(self,theta):
        
        
        assert np.all(theta > 0) and np.all(theta < 1)
        powr = (1+self.pars['couple_rts'])/(self.pars['couple_rts']+self.pars['crra_power'])
        tf = theta
        tm = 1-theta
        ces = (tf**powr + tm**powr)**(1/powr)
        umult = (self.pars['A']**(1-self.pars['crra_power']))*ces
        
        
        
        assert umult.shape == theta.shape
        
        return umult
    
    # this is multipliers for consumption relative to couple's level of consumption 
    def c_mult_just_c(self,theta):
        
        
        assert np.all(theta > 0) and np.all(theta < 1)
        powr = (1+self.pars['couple_rts'])/(self.pars['couple_rts']+self.pars['crra_power'])
        irho = 1/(1+self.pars['couple_rts'])
        irs  = 1/(self.pars['couple_rts']+self.pars['crra_power'])
        tf = theta
        tm = 1-theta
        bottom = (tf**(powr) + tm**(powr))**irho 
        
        kf = self.pars['A']*(tf**(irs))/bottom
        km = self.pars['A']*(tm**(irs))/bottom
        
        assert kf.shape == theta.shape
        assert km.shape == theta.shape
        
        return kf, km
    
    
    def u_mult(self,theta,ind_l):
        
        # get raw umult
        
        umult_c = self.u_mult_just_c(theta)   
        # combine it with ls        
        u1l = self.ls_u1[ind_l]   
        sigma = self.pars['crra_power'] 
        oos = 1/sigma
        umult = (u1l**oos + umult_c**oos)**sigma
        
        
        #assert umult.shape == theta.shape
        
        return umult
        
    def c_mult(self,theta,ind_l):
        u1l = self.ls_u1[ind_l]
        umult_c = self.u_mult_just_c(theta)
        
        
        sigma = self.pars['crra_power'] 
        
        oos = 1/sigma
        
        k_x = u1l**oos / (u1l**oos + umult_c**oos)
        k_c = 1-k_x
        kf, km = self.c_mult_just_c(theta)
        
        k_cf = k_c*kf
        k_cm = k_c*km
        
        return k_cf, k_cm, k_x
        
        
        
    def u(self,c):
        return u_aux(c,self.pars['crra_power'])#(c**(1-self.pars['crra_power']))/(1-self.pars['crra_power'])
    
    def u_part(self,c,theta,ind_l,psi,ushift): # this returns utility of each partner out of some c
    
        k_cf, k_cm, k_x = self.c_mult(theta,ind_l)
        if np.any(k_x < 1e-5):
            assert np.all(k_x<1e-5)  # ! this makes sure that there are no theta
            # such that x = 0 with one theta and x > 0 with another
            u_both = self.ls_u0[ind_l] + psi + ushift
        else:
            u_both = self.ls_u1[ind_l]*self.u(k_x*c) + self.ls_u0[ind_l] + psi + ushift
        uf, um = self.u(k_cf*c) + u_both, self.u(k_cm*c) + u_both
        
        return uf, um
    
    def u_couple(self,c,theta,ind_l,psi,ushift): # this returns utility of each partner out of some c
        
        umult = self.u_mult(theta,ind_l)
        
            
        return umult*self.u(c) + self.ls_u0[ind_l] + psi + ushift
    
    
    def vm_last(self,s,zm,zf,psi,theta,ushift):
        # this is the value function for couple that has savings s,
        # Z = (zm,zf,psi) and bargaining power theta after all decisions are made
        
        
        #Get utility for different FLS
        #TODO check indeces here
        u_couple_g=np.zeros((self.na,self.pars['nexo_t'][-1],len(self.thetagrid),len(self.ls_levels)),dtype=self.dtype)
        income_g=np.zeros((self.na,self.pars['nexo_t'][-1],len(self.thetagrid),len(self.ls_levels)),dtype=self.dtype)
        u_couple_g=np.zeros((self.na,self.pars['nexo_t'][-1],len(self.thetagrid),len(self.ls_levels)),dtype=self.dtype)
        u_f_g=np.zeros((self.na,self.pars['nexo_t'][-1],len(self.thetagrid),len(self.ls_levels)),dtype=self.dtype)
        u_m_g=np.zeros((self.na,self.pars['nexo_t'][-1],len(self.thetagrid),len(self.ls_levels)),dtype=self.dtype)
        
        ftrend = self.pars['f_wage_trend'][-1]
        mtrend = self.pars['m_wage_trend'][-1]
        
        for l in range(len(self.ls_levels)):
           
            income_g[...,l] = self.pars['R_t'][-1]*s + np.exp(zm+mtrend) +  np.exp(zf+ftrend)*self.ls_levels[l]
            u_couple_g[...,l] = self.u_couple(income_g[...,l],theta,l,psi,ushift)
            u_f_g[...,l], u_m_g[...,l] = self.u_part(income_g[...,l],theta,l,psi,ushift)
        #Get optimal FLS
        lsi=np.argmax(u_couple_g,axis=3)[...,None]
        V  = np.take_along_axis(u_couple_g,lsi,axis=3).squeeze(axis=3)#[:,:,:,0]
        VF = np.take_along_axis(u_f_g,lsi,axis=3).squeeze(axis=3)
        VM = np.take_along_axis(u_m_g,lsi,axis=3).squeeze(axis=3)
        income=np.take_along_axis(income_g,lsi,axis=3).squeeze(axis=3)#[:,:,:,0]        
        lsout = self.ls_levels[lsi]
        
        
        return V.astype(self.dtype), VF.astype(self.dtype), VM.astype(self.dtype), income.astype(self.dtype), np.zeros_like(income.astype(self.dtype)), lsout.astype(self.dtype), u_couple_g.astype(self.dtype)
        

    def vm_last_grid(self,ushift):
        # this returns value of vm on the grid corresponding to vm
        s_in = self.agrid_c[:,None,None]
        zm_in = self.exogrid.all_t[-1][:,1][None,:,None]
        zf_in = self.exogrid.all_t[-1][:,0][None,:,None]
        psi_in = self.exogrid.all_t[-1][:,2][None,:,None]
        theta_in = self.thetagrid[None,None,:]
                
        return self.vm_last(s_in,zm_in,zf_in,psi_in,theta_in,ushift)
        
    
    

    def vs_last(self,s,z_plus_trend,ushift,return_cs=False):  
        # generic last period utility for single agent
        income = self.pars['R_t'][-1]*s+np.exp(z_plus_trend) 
        if return_cs:
            return self.u(income).astype(self.dtype) + ushift, income.astype(self.dtype), np.zeros_like(income.astype(self.dtype))
        else:
            return self.u(income)
    
    def vs_last_grid(self,female,ushift,return_cs=False):
        # this returns value of vs on the grid corresponding to vs
        s_in = self.agrid_s[:,None]
        z_in = self.exogrid.zf_t[-1][None,:] if female else self.exogrid.zm_t[-1][None,:]
        trend = self.pars['f_wage_trend'][-1] if female else self.pars['m_wage_trend'][-1]        
        return self.vs_last(s_in,z_in+trend,ushift,return_cs)
        
    
    
    
    
    
    

#from numba import jit
#@jit(nopython=True)
def u_aux(c,sigma):
    # this is pretty important not to have (c^sigma - 1) here as it is hard to 
    # keep everywhere and occasionally this generates nasty mistakes
    if sigma!=1:
        return (c**(1-sigma))/(1-sigma)
    else:
        return np.log(c)

    


class DivorceCosts(object):
    # this is something that regulates divorce costs
    # it aims to be fully flexible
    def __init__(self, 
                 unilateral_divorce=True, # whether to allow for unilateral divorce
                 assets_kept = 1.0, # how many assets of couple are splited (the rest disappears)
                 u_lost_m=0.0,u_lost_f=0.0, # pure utility losses b/c of divorce
                 money_lost_m=0.0,money_lost_f=0.0, # pure money (asset) losses b/c of divorce
                 money_lost_m_ez=0.0,money_lost_f_ez=0.0, # money losses proportional to exp(z) b/c of divorce
                 eq_split=0.0 #The more of less equal way assets are split within divorce
                 ): # 
        
        self.unilateral_divorce = unilateral_divorce # w
        self.assets_kept = assets_kept
        self.u_lost_m = u_lost_m
        self.u_lost_f = u_lost_f
        self.money_lost_m = money_lost_m
        self.money_lost_f = money_lost_f
        self.money_lost_m_ez = money_lost_m_ez
        self.money_lost_f_ez = money_lost_f_ez
        self.eq_split = eq_split
        
    def shares_if_split(self,income_share_f):
        
        
        shf=(0.5*self.eq_split + income_share_f*(1-self.eq_split))
        share_f = self.assets_kept*shf - self.money_lost_f
        share_m = self.assets_kept*(1-shf) - self.money_lost_m
        
        return share_f, share_m
        
        

def build_s_grid(agrid,n_between,da_min,da_max):
    sgrid = np.array([0.0],np.float64)
    for j in range(agrid.size-1):
        step = (agrid[j+1] - agrid[j])/n_between
        if step >= da_min and step <= da_max:
            s_add = np.linspace(agrid[j],agrid[j+1],n_between)[:-1]
        elif step < da_min:
            s_add = np.arange(agrid[j],agrid[j+1],da_min)
        elif step > da_max:
            s_add = np.arange(agrid[j],agrid[j+1],da_max)
        sgrid = np.concatenate((sgrid,s_add))
    
    sgrid = np.concatenate((sgrid,np.array([agrid[-1]])))
            
    if sgrid[0] == sgrid[1]: 
        sgrid = sgrid[1:]
        
    return sgrid
