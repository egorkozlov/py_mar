#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This contains things relevant for setting up the model
"""

import numpy as np

from rw_approximations import rouw_nonst, normcdf_tr,tauchen_nonst
from mc_tools import combine_matrices_two_lists, int_prob,cut_matrix
from scipy.stats import norm
from scipy import optimize
from collections import namedtuple
from gridvec import VecOnGrid
import pickle
from scipy import sparse



class ModelSetup(object):
    def __init__(self,nogrid=False,divorce_costs='Default',separation_costs='Default',**kwargs): 
        p = dict()       
        period_year=1#this can be 1,2,3 or 6
        transform=2#this tells how many periods to pull together for duration moments
        T = int(62/period_year)#int(62/period_year)
        Tret = int(42/period_year)#int(42/period_year) # first period when the agent is retired
        Tbef=int(0/period_year)
        Tren  = int(42/period_year)#int(42/period_year)#int(42/period_year) # period starting which people do not renegotiate/divroce
        Tmeet = int(42/period_year)#int(42/period_year)#int(42/period_year) # period starting which you do not meet anyone
        p['py']=period_year
        p['ty']=transform
        p['T'] = T
        p['Tret'] = Tret
        p['Tren'] = Tren
        p['Tbef'] = Tbef
        p['sig_zf_0']  = 0.5449176#0.4096**(0.5)
        p['sig_zf']    = .0277534**(0.5)#0.0399528**(0.5)
        p['sig_zm_0']  = .4435924#.405769**(0.5)
        p['sig_zm']    = .025014**(0.5)#0.0417483**(0.5)
        p['n_zf_t']      = [4]*Tret + [4]*(T-Tret)
        p['n_zm_t']      = [3]*Tret + [3]*(T-Tret)
        p['n_zf_correct']=1
        p['sigma_psi_mult'] = 0.28
        p['sigma_psi']   = 0.11
        p['n_psi_t']     = [15]*T
        p['R_t'] = [1.02**period_year]*T
        p['beta_t'] = [0.98**period_year]*T
        p['A'] = 1.0 # consumption in couple: c = (1/A)*[c_f^(1+rho) + c_m^(1+rho)]^(1/(1+rho))
        p['crra_power'] = 1.5
        p['couple_rts'] = 0.0 
        p['sig_partner_a'] = 0.1#0.05
        p['sig_partner_zf'] = 1.5#4.0#1.8#0.4 #This is crazy powerful for the diff in diff estimate
        p['sig_partner_zm'] = 0.95#0.5
        p['sig_partner_mult'] = 1.0
        p['dump_factor_z'] = 0.25#0.85
        p['dump_factor_a'] = 0.8#0.65
        p['mean_partner_z_female'] = 0.01#0.05#0.02
        p['mean_partner_z_male'] =  -0.03#0.01#-0.02
        p['mean_partner_a_female'] = 0.22#0.32
        p['mean_partner_a_male'] = -0.22#-0.32
        p['m_bargaining_weight'] = 0.5
        p['pmeet'] = 0.5
        
        p['z_drift'] = -0.09#-0.1
        p['rho_s']    =  0.0 #correlation in shocks
        
        p['wage_gap'] = 0.6
        p['wret'] = 0.6#0.5
        p['uls'] = 0.2
        p['pls'] = 1.0
        
        
        
        p['u_shift_mar'] = 0.0
        p['u_shift_coh'] = 0.0#-0.17625478002929690
        
         
        #Wages over time
        p['f_wage_trend'] = [0.0*(t>=Tret)+(t<Tret)*(-.74491918 +.04258303*(t-2) -.0016542*(t-2)**2+.00001775*(t-2)**3) for t in range(T)]
        p['f_wage_trend_single'] = [0.0*(t>=Tret)+(t<Tret)*(-.6805060 +.04629912*(t-2) -.00160467*(t-2)**2+.00001626*(t-2)**3) for t in range(T)]
        p['m_wage_trend'] = [0.0*(t>=Tret)+(t<Tret)*(-0.5620125  +0.06767721*t -0.00192571*t**2+ 0.00001573*t**3) for t in range(T)]
        p['m_wage_trend_single'] = [0.0*(t>=Tret)+(t<Tret)*(-.5960803  +.05829568*t -.00169143*t**2+ .00001446*t**3) for t in range(T)]
   
              
        p['f_wage_trend'] = [0.0*(t>=Tret+2)+(t<Tret+2)*(0.00+-.77138877  +.05915875*t -.00232914*t**2+ .00002484*t**3) for t in range(T)]
        p['f_wage_trend_single'] =  [0.0*(t>=Tret+2)+(t<Tret+2)*(0.00+-.67980802  +.04603417*t -.00158584*t**2+ .00001594*t**3) for t in range(T)]
        p['m_wage_trend'] = [0.0*(t>=Tret)+(t<Tret)*(-.434235  +.06016318*t -.00183131*t**2+ .00001573*t**3) for t in range(T)]
        p['m_wage_trend_single'] = [0.0*(t>=Tret)+(t<Tret)*(-.486139  +.05170349*t -.00160466*t**2+ .00001446*t**3) for t in range(T)]
           
  
        #Build trend for single assets
        p['ass_ratio']=[max(0,9.813475 -1.353955*t+ 0.0655805*t**2- 0.0007908*t**3)  for t in range(T)]
        p['ass_f']=p['ass_ratio']*np.exp(p['f_wage_trend_single'])
        p['ass_m']=p['ass_ratio']*np.exp(p['m_wage_trend_single'])
        
        
        p['util_lam'] = 0.189#0.4
        p['util_alp_temp'] = 0.5
        p['util_xi'] = 1.07
        p['util_kap_temp']=0.206
        p['rprice_durables'] = 1.0#
        
        
        with open('assets.pkl', 'rb') as file:assets=pickle.load(file)
        p['av_a_m']=assets['av_a_m']
        p['av_a_f']=assets['av_a_f']
        p['totf']=assets['totf']
        p['totm']=assets['totm']
       
        

        
        for key, value in kwargs.items():
            assert (key in p), 'wrong name?'
            p[key] = value
            
            
        p['util_kap'] = (1-p['util_kap_temp'])/(p['util_kap_temp'])
        
        #Adjust aplpha do be like in greenoood
        p['util_alp']=p['util_alp_temp']*(p['util_kap_temp'])**((1-p['util_xi'])/p['util_lam'])
        
        #Adjust kappa and alpha to make sense of relative prices
        p['util_alp_m']=p['util_alp']*(1.0/(p['rprice_durables'])**(1.0-p['util_xi']))
        p['util_kap_m']=p['util_kap']*p['rprice_durables']**p['util_lam']
            
            
        # no replacements after this pint     
        p['sigma_psi_init'] = p['sigma_psi_mult']*p['sigma_psi']
        
        #Get the probability of meeting, adjusting for year-period
        p_meet=p['pmeet']
        for j in range(period_year-1):
            p_meet=p_meet+(1-p_meet)*p['pmeet']
            
            
        # timing here    
            
        p['pmeet_t'] = [p_meet]*Tmeet + [0.0]*(T-Tmeet)
        p['can divorce'] = [True]*Tren + [False]*(T-Tren)
        
        
        self.pars = p
        
        self.dtype = np.float64 # type for all floats
        
       
        
        # relevant for integration
        self.state_names = ['Female, single','Male, single','Couple, M', 'Couple, C']
        
        # female labor supply
        self.ls_levels = np.array([0.0,0.811],dtype=self.dtype)
        self.mlevel=0.811
        #self.ls_utilities = np.array([p['uls'],0.0],dtype=self.dtype)
        self.ls_pdown = np.array([p['pls'],0.0],dtype=self.dtype)
        self.nls = len(self.ls_levels)
        
        
        
        #Trim values
        self.trim=0.00000000001
        
        
        #Cost of Divorce
        if divorce_costs == 'Default':
            # by default the costs are set in the bottom
            self.div_costs = DivorceCosts(eq_split=1.0,assets_kept=1.0)
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
            self.sep_costs = DivorceCosts(eq_split=0.0,assets_kept=1.0)
        else:
            if isinstance(separation_costs,dict):
                # you can feed in arguments to DivorceCosts
                self.sep_costs = DivorceCosts(**separation_costs)
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
            
            exogrid['zf_t'],  exogrid['zf_t_mat'],zft,zftmat,exogrid['zm_t'],  exogrid['zm_t_mat']=dict(),dict(),dict(),dict(),dict(),dict()
            zft,       zftmat                    = rouw_nonst(p['T'],p['sig_zf']*period_year**0.5,p['sig_zf_0'],p['n_zf_t'][0]-p['n_zf_correct'])
            exogrid['zm_t'],  exogrid['zm_t_mat']= rouw_nonst(p['T'],p['sig_zm']*period_year**0.5,p['sig_zm_0'],p['n_zm_t'][0])
            #Embody the grid for women in a bigger one
            if p['n_zf_correct']>0:

                exogrid['zf_t']=list()
                exogrid['zf_t_mat']=list()
                for t in range(p['T']):
                    
                    
                    #Extend grid
                    h=zft[t][1]-zft[t][0]
                    # dist1=zft[t][0]-h
                    # dist0=zft[t][0]-p['n_zf_correct']*h
                    dist2=zft[t][0]-h
                    dist1=zft[t][0]-(p['n_zf_correct']-1)*h
                    dist0=zft[t][0]-p['n_zf_correct']*h
                    
                    #Copy transition matrix
                    #exogrid['zf_t']=exogrid['zf_t']+[np.concatenate((np.array([dist0,dist1,dist2]),zft[t]))]
                    #exogrid['zf_t']=exogrid['zf_t']+[np.concatenate((np.array([dist0,dist1]),zft[t]))]
                    exogrid['zf_t']=exogrid['zf_t']+[np.concatenate((np.array([dist0]),zft[t]))]
                    exogrid['zf_t_mat']=exogrid['zf_t_mat']+[np.zeros((p['n_zf_t'][t],p['n_zf_t'][t]))]
                    exogrid['zf_t_mat'][t][p['n_zf_correct']:,p['n_zf_correct']:]=zftmat[t]
                    
                    #Shift transition matrix to fill values
                    if t<p['T']-1:
                        
                        exogrid['zf_t_mat'][t][0,:-p['n_zf_correct']]=zftmat[t][0,:]
                        #exogrid['zf_t_mat'][t][1,:-p['n_zf_correct']]=zftmat[t][1,:]
                        #exogrid['zf_t_mat'][t][2,:-p['n_zf_correct']]=zftmat[t][2,:]
                       
                            
                    else:
                        exogrid['zf_t_mat'][t]=None
                       
                    
            else:    

                exogrid['zf_t']=zft
                exogrid['zf_t_mat']=zftmat
                       
                    


                    
                    
            #Drift the grids
            for t in range(Tret):
                exogrid['zf_t'][t]=exogrid['zf_t'][t]
            for t in range(Tret-2):
                exogrid['zm_t'][t]=exogrid['zm_t'][t]
                    
                    
            
            ################################
            #First mimic US pension system
            ###############################
            

            
            #function to compute pension
            def pens(value):
                
                #Get median income before retirement using men model income in Tret-1
                #+ ratio of men and female labor income for the rest
                yret=(1.73377+(.8427056/1.224638)*1.73377* 0.3246206)/(1+0.3246206)
                thresh1=0.38*yret
                thresh2=1.59*yret
                
                inc1=np.minimum(value,thresh1)
                inc2=np.maximum(np.minimum(value-inc1,thresh2-inc1),0)
                inc3=np.maximum(value-thresh2,0)
                
                return inc1*0.9+inc2*0.32+inc3*0.15
                
              
            
            # for t in range(Tret,T):
            #     exogrid['zf_t'][t] = np.array([np.log(p['wret'])])
            #     exogrid['zm_t'][t] = np.array([np.log(p['wret'])])
            #     exogrid['zf_t_mat'][t] = np.atleast_2d(1.0)
            #     exogrid['zm_t_mat'][t] = np.atleast_2d(1.0)
                
                
            # # fix transition from non-retired to retired    
            # exogrid['zf_t_mat'][Tret-1] = np.ones((p['n_zf_t'][Tret-1],1))
            # exogrid['zm_t_mat'][Tret-1] = np.ones((p['n_zm_t'][Tret-1],1))
            
            # #Tax system as in Wu and Kruger
            # for t in range(0,Tret):
            #     exogrid['zf_t'][t] = exogrid['zf_t'][t]#*(1-0.1327)+np.log(1-0.1575)
            #     exogrid['zm_t'][t] = exogrid['zm_t'][t]#*(1-0.1327)+np.log(1-0.1575)  
            
            #Comment out the following if you dont want retirment based on income
            for t in range(Tret+2,T):
                exogrid['zf_t'][t] = np.log(pens(np.exp(p['f_wage_trend'][Tret+1]+exogrid['zf_t'][Tret+1])))#np.array([np.log(p['wret'])])
                             
                exogrid['zf_t_mat'][t] = np.diag(np.ones(len(exogrid['zf_t'][t])))#p.atleast_2d(1.0)
                
                
            for t in range(Tret,T):    
                exogrid['zm_t'][t] = np.log(pens(np.exp(p['m_wage_trend'][Tret-1]+exogrid['zm_t'][Tret-1])))  
                exogrid['zm_t_mat'][t] = np.diag(np.ones(len(exogrid['zm_t'][t])))
        
                
            
            # fix transition from non-retired to retired    
            exogrid['zf_t_mat'][Tret+1] = np.diag(np.ones(len(exogrid['zf_t'][Tret+1])))
            exogrid['zm_t_mat'][Tret-1] = np.diag(np.ones(len(exogrid['zm_t'][Tret-1])))


            
            exogrid['psi_t'], exogrid['psi_t_mat'] = tauchen_nonst(p['T'],p['sigma_psi']*period_year**0.5,p['sigma_psi_init'],p['n_psi_t'][0])
            exogrid['psi_t_mat'][Tret-1] = np.diag(np.ones(len(exogrid['psi_t_mat'][Tret-1])))
            exogrid['psi_t_mat'][Tret] = np.diag(np.ones(len(exogrid['psi_t_mat'][Tret-1])))
            exogrid['psi_t_mat'][Tret+1] = np.diag(np.ones(len(exogrid['psi_t_mat'][Tret-1])))
            
   
            # #Here I impose no change in psi from retirement till the end of time 
            # for t in range(Tret,T-1):
               
            #     exogrid['psi_t'][t] = exogrid['psi_t'][Tret-1]#np.array([np.log(p['wret'])])             
            #     exogrid['psi_t_mat'][t] = np.diag(np.ones(len(exogrid['psi_t'][t])))#p.atleast_2d(1.0)

            zfzm, zfzmmat = combine_matrices_two_lists(exogrid['zf_t'], exogrid['zm_t'], exogrid['zf_t_mat'], exogrid['zm_t_mat'],trim=self.trim)

            
            #Create a new bad version of transition matrix p(zf_t)
            
            
            zf_bad = [tauchen_drift(exogrid['zf_t'][t].copy(), exogrid['zf_t'][t+1].copy(), 
                                                1.0, p['sig_zf'], p['z_drift'], exogrid['zf_t_mat'][t])
                                    for t in range(self.pars['Tret']-1) ]
            
            #Account for retirement here
            zf_bad = zf_bad+[exogrid['zf_t_mat'][t] for t in range(self.pars['Tret']-1,self.pars['T']-1)]+ [None]
            
            zf_t_mat_down = zf_bad
            
            zfzm2, zfzmmat2 = combine_matrices_two_lists(exogrid['zf_t'], exogrid['zm_t'], zf_t_mat_down, exogrid['zm_t_mat'],trim=self.trim)
            
            if p['rho_s']>0:
                for t in range(p['Tret']-1):
                    for j in range(p['n_zm_t'][t]):
                        for ym in range(p['n_zm_t'][t]):
                        
                            
                            rhom=(1.0-p['rho_s']**2)**0.5
                            prec=exogrid['zm_t'][t][j] if t>0 else 0.0
                            drif=p['rho_s']*p['sig_zf']/p['sig_zm']*(exogrid['zm_t'][t+1][ym]-prec)
                            mat1=tauchen_drift(exogrid['zf_t'][t].copy(), exogrid['zf_t'][t+1].copy(), 1.0, rhom*p['sig_zf'], drif, exogrid['zf_t_mat'][t])
                            mat2=tauchen_drift(exogrid['zf_t'][t].copy(), exogrid['zf_t'][t+1].copy(), 1.0, rhom*p['sig_zf'], drif+p['z_drift'], exogrid['zf_t_mat'][t])
                            for i in range(p['n_zf_t'][t]): 
                        
                                #Modify the grid for women
                                exogrid['zf_t_mat2'][t][i,:]= mat1[i,:]

                                exogrid['zf_t_mat2d'][t][i,:]=mat2[i,:]
                                
                                ##Update the big Matrix
                                for yf in range(p['n_zf_t'][t]):
                                
                                    
                                    zfzmmat[t][i*(p['n_zm_t'][t]-1)+j+i,yf*(p['n_zm_t'][t]-1)+ym+yf]=exogrid['zf_t_mat2'][t][i,yf]*exogrid['zm_t_mat'][t][j,ym]
                                    zfzmmat2[t][i*(p['n_zm_t'][t]-1)+j+i,yf*(p['n_zm_t'][t]-1)+ym+yf]=exogrid['zf_t_mat2d'][t][i,yf]*exogrid['zm_t_mat'][t][j,ym]
           
            
            #Adjust retirement as in Heatcote et al.
            for t in range(p['Tret'],p['Tret']+2):
                for j in range(len(zfzm[t])):
                    pref=np.exp(zfzm[t][j][0])+np.exp(zfzm[t][j][1])
                    zfzm[t][j][0]=-20.0
                    zfzm[t][j][1]=np.log(pref)
                    pref=np.exp(zfzm2[t][j][0])+np.exp(zfzm2[t][j][1])
                    zfzm2[t][j][0]=-20.0
                    zfzm2[t][j][1]=np.log(pref)
                    
                    
            for t in range(p['Tret']+2,p['T']):
                for j in range(len(zfzm[t])):
                    pref=max(np.exp(zfzm[t][j][0])+np.exp(zfzm[t][j][1]),1.5*max(np.exp(zfzm[t][j][0]),np.exp(zfzm[t][j][1])))
                    zfzm[t][j][0]=-20.0
                    zfzm[t][j][1]=np.log(pref)
                    pref=max(np.exp(zfzm2[t][j][0])+np.exp(zfzm2[t][j][1]),1.5*max(np.exp(zfzm2[t][j][0]),np.exp(zfzm2[t][j][1])))
                    zfzm2[t][j][0]=-20.0
                    zfzm2[t][j][1]=np.log(pref)
            
            
            #Put everything together
            all_t, all_t_mat = combine_matrices_two_lists(zfzm,exogrid['psi_t'],zfzmmat,exogrid['psi_t_mat'],trim=self.trim)
            all_t_mat_sparse_T = [sparse.csc_matrix(D.T) if D is not None else None for D in all_t_mat]
            
            all_t_down, all_t_mat_down = combine_matrices_two_lists(zfzm2,exogrid['psi_t'],zfzmmat2,exogrid['psi_t_mat'],trim=self.trim)
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
        self.na = 40#40
        self.scala=1.0
        self.amin = 0
        self.amax = 60*self.scala
        self.amax1 = 100*self.scala
        self.agrid_c = np.linspace(self.amin,self.amax,self.na,dtype=self.dtype)
        tune=10#30.5
        self.agrid_c = np.geomspace(self.amin+tune,self.amax+tune,num=self.na)-tune
        self.agrid_c[-1]=self.amax1
        self.agrid_c[-2]=80*self.scala
        # this builds finer grid for potential savings
        s_between = 7 # default numer of points between poitns on agrid
        s_da_min = 0.001*self.scala # minimal step (does not create more points)
        s_da_max = 5.0*self.scala # maximal step (creates more if not enough)
        
        self.sgrid_c = build_s_grid(self.agrid_c,s_between,s_da_min,s_da_max)
        self.vsgrid_c = VecOnGrid(self.agrid_c,self.sgrid_c)
        
        
         
        #Grid Single
        scale = 1.1
        self.amin_s = 0
        self.amax_s = self.amax/scale
        self.agrid_s = np.linspace(self.amin_s,self.amax_s,self.na,dtype=self.dtype)
        #self.agrid_s[self.na-1]=18#180
        tune_s=2.5
        self.agrid_s = self.agrid_c/2#np.geomspace(self.amin_s+tune_s,self.amax_s+tune_s,num=self.na)-tune_s
        #self.agrid_s[-1]=self.amax1/scale
        #self.agrid_c[-2]=120/scale
        self.sgrid_s = build_s_grid(self.agrid_s,s_between,s_da_min,s_da_max)
        self.vsgrid_s = VecOnGrid(self.agrid_s,self.sgrid_s)
        
        # grid for theta
        self.ntheta = 13
        self.thetamin = 0.1
        self.thetamax = 0.9
        self.thetagrid = np.linspace(self.thetamin,self.thetamax,self.ntheta,dtype=self.dtype)
        
        #Grid for the share in assets
        self.ashare = np.linspace(0.05,0.95,3,dtype=self.dtype)#self.ashare = np.linspace(0.15,0.85,3,dtype=self.dtype)
        
        
        
        
        
        
        # construct finer grid for bargaining
        ntheta_fine = 1*self.ntheta # actual number may be a bit bigger
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
        
        #Get indexes from fine back to coarse thetagrid
        cg=VecOnGrid(self.thetagrid,self.thetagrid_fine)
        index_t=cg.i
        index_t1=index_t+1
        wherep=(cg.wthis<0.5)
        self.igridcoarse=index_t
        self.igridcoarse[wherep]=index_t1[wherep]


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
        
        
        # building m grid
        ezfmin = min([np.min(self.ls_levels[-1]*np.exp(g+t)) for g,t in zip(exogrid['zf_t'],p['f_wage_trend'])])
        ezmmin = min([np.min(self.mlevel*np.exp(g+t)) for g,t in zip(exogrid['zm_t'],p['m_wage_trend'])])
        ezfmax = max([np.max(self.ls_levels[-1]*np.exp(g+t)) for g,t in zip(exogrid['zf_t'],p['f_wage_trend'])])
        ezmmax = max([np.max(self.mlevel*np.exp(g+t)) for g,t in zip(exogrid['zm_t'],p['m_wage_trend'])])
        
        
        
        self.money_min = 0.95*min(ezmmin,ezfmin) # cause FLS can be up to 0
        #self.mgrid = ezmmin + self.sgrid_c # this can be changed later
        mmin = self.money_min
        mmax = ezfmax + ezmmax + np.max(self.pars['R_t'])*self.amax1
        mint = (ezfmax + ezmmax) # poin where more dense grid begins
        
        ndense = 1200
        nm = 3000
        
        gsparse = np.linspace(mint,mmax,nm-ndense)
        gdense = np.linspace(mmin,mint,ndense+1) # +1 as there is a common pt
        
        
       
        
        self.mgrid = np.zeros(nm,dtype=self.dtype)
        self.mgrid[ndense:] = gsparse
        self.mgrid[:(ndense+1)] = gdense
        #self.mgrid = np.geomspace(mmin+tune,mmax+tune,num=nm)-tune#build_s_grid(self.sgrid_c,10,s_da_min*0.1,s_da_max*0.01)+mmin#
        
        assert np.all(np.diff(self.mgrid)>0)
        
        self.u_precompute()
        
        
    def mar_mats_assets(self,npoints=12,abar=0.1):
        # for each grid point on single's grid it returns npoints positions
        # on (potential) couple's grid's and assets of potential partner 
        # (that can be off grid) and correpsonding probabilities. 
        
        self.prob_a_mat = dict()
        self.i_a_mat = dict()
        
        na = self.agrid_s.size
        
        agrid_s = self.agrid_s
        agrid_c = self.agrid_c
        
        s_a_partner = self.pars['sig_partner_a']
        
        for female in [True,False]:
            prob_a_mat = np.zeros((self.pars['T'],na,npoints),dtype=self.dtype)
            i_a_mat = np.zeros((self.pars['T'],na,npoints),dtype=np.int16)
            mena=-0.15 if female else 0.15
            
            
            for ia, a in enumerate(agrid_s):
                for t in range(self.pars['T']):
                    lagrid_t = np.zeros_like(agrid_c)
                    
                    i_neg = (agrid_c <= max(abar,a-mena) - 1e-6)
                    
                    # if a is zero this works a bit weird but does the job
                    
                    lagrid_t[~i_neg] = np.log(2e-6 + (agrid_c[~i_neg] - max((a-mena),0))/max(abar,a))#agrid_c[~i_neg]#
                    
                    
                    lmin = lagrid_t[~i_neg].min()
                    # just fill with very negative values so this is never chosen
                    lagrid_t[i_neg] = lmin - s_a_partner*10 - \
                        s_a_partner*np.flip(np.arange(i_neg.sum())) 
                    
                    # TODO: this needs to be checked
                    if female:
                        me=self.pars['av_a_m'][t]
                        mean=self.pars['mean_partner_a_female']#a*1.6+me*0.4#np.log(2e-6 + ( me)/max(abar,a))#
                        st=2.5*max(np.std((2e-6 + (self.pars['totm'][t]))),0.01)#max(np.std(np.log(2e-6 + (self.pars['totm'][t])/max(abar,a))),0.001)#s_a_partner
                    else:
                        me=self.pars['av_a_f'][t]
                        mean=self.pars['mean_partner_a_male']#(a*1.6*0.95+me*0.4)#np.log(2e-6 + (me )/max(abar,a))#
                        st=max(np.std((2e-6 + (self.pars['totf'][t]))),0.01)#max(np.std(np.log(2e-6 + (self.pars['totf'][t])/max(abar,a))),0.001)#s_a_partner
                        
                    #p_a = int_prob(lagrid_t,mu=mean,sig=st,n_points=npoints)
                    p_a = int_prob(lagrid_t,mu=mean,sig=st**0.5,n_points=npoints)
                    
                    
                    p_a  = int_prob(lagrid_t, mu=self.pars['dump_factor_a']*mean
                                      ,sig=(1-self.pars['dump_factor_a'])**
                                      0.5*s_a_partner*self.pars['sig_partner_mult'],n_points=npoints)
                     
                    i_pa = (-p_a).argsort()[:npoints] # this is more robust then nonzero
                    p_pa = p_a[i_pa]
                    prob_a_mat[t,ia,:] = p_pa
                    i_a_mat[t,ia,:] = i_pa
            
            
            self.prob_a_mat[female] = prob_a_mat
            self.i_a_mat[female] = i_a_mat
            

        
    
    def mar_mats_iexo(self,t,female=True):
        # TODO: check timing
        # this returns transition matrix for single agents into possible couples
        # rows are single's states
        # columnts are couple's states
        # you have to transpose it if you want to use it for integration
        setup = self
        trim_lvl=setup.trim
        
        nexo = setup.pars['nexo_t'][t]
        sigma_psi_init = setup.pars['sigma_psi_init']
        #sig_z_partner = setup.pars['sig_partner_z']
        psi_couple = setup.exogrid.psi_t[t+1]
        
        
        if female:
            nz_single = setup.exogrid.zf_t[t].shape[0]
            p_mat = np.empty((nexo,nz_single))
            z_own = setup.exogrid.zf_t[t]
            n_zown = z_own.shape[0]
            z_partner = setup.exogrid.zm_t[t]
            zmat_own = setup.exogrid.zf_t_mat[t]
            trend=setup.pars['m_wage_trend_single'][t]
            mean=setup.pars['mean_partner_z_female']-setup.pars['m_wage_trend'][t]+setup.pars['m_wage_trend_single'][t]
            sig_z_partner=(setup.pars['sig_zm_0']**2+(t+1)*setup.pars['sig_partner_zm']*setup.pars['sig_zm']**2)**0.5
        else:
            nz_single = setup.exogrid.zm_t[t].shape[0]
            p_mat = np.empty((nexo,nz_single))
            z_own = setup.exogrid.zm_t[t]
            n_zown = z_own.shape[0]
            z_partner = setup.exogrid.zf_t[t]
            zmat_own = setup.exogrid.zm_t_mat[t]    
            trend=setup.pars['f_wage_trend_single'][t]
            mean=setup.pars['mean_partner_z_male']-setup.pars['f_wage_trend'][t]+setup.pars['f_wage_trend_single'][t]
            sig_z_partner=(setup.pars['sig_zf_0']**2+(t+1)*setup.pars['sig_partner_zf']*setup.pars['sig_zf']**2)**0.5
            
        def ind_conv(a,b,c): return setup.all_indices(t,(a,b,c))[0]
        
        
        for iz in range(n_zown):
            p_psi = int_prob(psi_couple,mu=0,sig=sigma_psi_init)
            if female:
                p_zm  = int_prob(z_partner, mu=setup.pars['dump_factor_z']*z_own[iz]+
                                  mean+setup.pars['mean_partner_z_female'],sig=(1-setup.pars['dump_factor_z'])**
                                  0.5*sig_z_partner*setup.pars['sig_partner_mult'])
                p_zf  = zmat_own[iz,:]
            else:
                p_zf  = int_prob(z_partner, mu=setup.pars['dump_factor_z']*z_own[iz]+ 
                                 mean+setup.pars['mean_partner_z_male'],sig=(1-setup.pars['dump_factor_z'])**
                                 0.5*sig_z_partner*setup.pars['sig_partner_mult'])
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
        
        

        
        self.matches = dict()
        
        for female in [True,False]:
            desc = 'Female, single' if female else 'Male, single'
            

            
            pmats = self.part_mats[desc] 
            
            
            match_matrix = list()
            
            
            for t in range(self.pars['T']-1):
                pmat_iexo = pmats[t] # nz X nexo
                # note that here we do not use transpose
                
                pmat_a = self.prob_a_mat[female][t,:,:]
                imat_a = self.i_a_mat[female][t,:,:]
            
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
    
    def u_mult(self,theta):
        assert np.all(theta > 0) and np.all(theta < 1)
        powr = (1+self.pars['couple_rts'])/(self.pars['couple_rts']+self.pars['crra_power'])
        tf = theta
        tm = 1-theta
        ces = (tf**powr + tm**powr)**(1/powr)
        umult = (self.pars['A']**(1-self.pars['crra_power']))*ces
        
        
        
        assert umult.shape == theta.shape
        
        return umult
    
    
    def c_mult(self,theta):
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
    
    def u(self,c):
        return u_aux(c,self.pars['crra_power'])#(c**(1-self.pars['crra_power']))/(1-self.pars['crra_power'])
    
    
    def u_pub(self,x,l,mt=0.0):
        alp = self.pars['util_alp_m']
        xi = self.pars['util_xi']
        lam = self.pars['util_lam']
        kap = self.pars['util_kap_m']        
        return alp*(x**lam + kap*(1+mt-l)**lam)**((1-xi)/lam)/(1-xi)
    
    
    def u_part(self,c,x,il,theta,ushift,psi): # this returns utility of each partner out of some c
        kf, km = self.c_mult(theta)   
        l = self.ls_levels[il]
        upub = self.u_pub(x,l,mt=1.0-self.mlevel) + ushift + psi
        return self.u(kf*c) + upub, self.u(km*c) + upub
    
    def u_couple(self,c,x,il,theta,ushift,psi): # this returns utility of each partner out of some c
        umult = self.u_mult(theta) 
        l = self.ls_levels[il]
        return umult*self.u(c) + self.u_pub(x,l,mt=1.0-self.mlevel) + ushift + psi
    
    def u_single_pub(self,c,x,l):
        return self.u(c) + self.u_pub(x,l)
    
        
        
    
    def u_precompute(self):
        from intratemporal import int_sol
        sig = self.pars['crra_power']
        alp = self.pars['util_alp_m']
        xi = self.pars['util_xi']
        lam = self.pars['util_lam']
        kap = self.pars['util_kap_m']
        
        nm = self.mgrid.size
        ntheta = self.ntheta
        nl = self.nls
        
        uout = np.empty((nm,ntheta,nl),dtype=self.dtype)
        xout = np.empty((nm,ntheta,nl),dtype=self.dtype)
        cout = np.empty((nm,ntheta,nl),dtype=self.dtype)
        
        for il in range(nl):
            for itheta in range(ntheta):
                A = self.u_mult(self.thetagrid[itheta])
                ls = self.ls_levels[il]
                x, c, u = int_sol(self.mgrid,A=A,alp=alp,sig=sig,xi=xi,lam=lam,kap=kap,lbr=ls,mt=1.0-self.mlevel)
                uout[:,itheta,il] = u
                xout[:,itheta,il] = x
                cout[:,itheta,il] = c
                
                
        self.ucouple_precomputed_u = uout
        self.ucouple_precomputed_x = xout
        self.ucouple_precomputed_c = cout
                
        
        # singles have just one level of labor supply (work all the time)
        
        xout, cout, uout = int_sol(self.mgrid,A=1,alp=alp,sig=sig,xi=xi,lam=lam,kap=kap,lbr=self.ls_levels[-1])#self.ls_levels[-1]
        self.usinglef_precomputed_u = uout
        self.usinglef_precomputed_x = xout
        self.usinglef_precomputed_c = cout
        xout, cout, uout = int_sol(self.mgrid,A=1,alp=alp,sig=sig,xi=xi,lam=lam,kap=kap,lbr=self.mlevel)
        self.usinglem_precomputed_u = uout
        self.usinglem_precomputed_x = xout
        self.usinglem_precomputed_c = cout
    

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
                 u_lost_m=0.0,u_lost_f=0.0, # pure utility losses b/c of divorceagri
                 money_lost_m=0.0,money_lost_f=0.0, # pure money (asset) losses b/c of divorce
                 money_lost_m_ez=0.0,money_lost_f_ez=0.0, # money losses proportional to exp(z) b/c of divorce
                 eq_split=1.0 #The more of less equal way assets are split within divorce
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
    
    
    def shares_if_split_theta(self,setup,theta,dc):
        
        #First build the title based sharing rule
        sharef=theta#setup.c_mult(theta)[0]
        shf=(0.5*dc.eq_split + sharef*(1-dc.eq_split))
        share_f = dc.assets_kept*shf
        share_m=dc.assets_kept*(1-shf)
        
        return share_f,share_m
    

       
        
def tauchen_drift(z_now,z_next,rho,sigma,mu,mat):
    z_now = np.atleast_1d(z_now)
    z_next = np.atleast_1d(z_next)
    if z_next.size == 1:
        return np.ones((z_now.size,1),dtype=z_now.dtype)
    
    d = np.diff(z_next)
    assert np.ptp(d) < 1e-5, 'Step size should be fixed'
    
    h_half = d[0]/2
    
    Pi = np.zeros((z_now.size,z_next.size),dtype=z_now.dtype)
    Pii = np.zeros((z_now.size,z_next.size),dtype=z_now.dtype)
    
    ez = rho*z_now + mu
    
    
    def f(x):
        
        pi=int_prob(z_next,mu=x,sig=sigma)
        return np.exp(ez[j])/np.exp(np.sum(z_next*pi))-1.0

    for j in range(z_next.size):
        Pi[j,:]=int_prob(z_next,mu=ez[j],sig=sigma)
        if (abs(ez[j]-np.sum(z_next*Pi[j,:]))>0.001):
            
            if (f(ez[j]-1.0)>0 and f(ez[j]+1.0)<0):
                sol = optimize.root_scalar(f, x0=ez[j],bracket=[ez[j]-1.0, ez[j]+1.0], maxiter=200,xtol=0.0001,method='bisect')
                mu1=sol.root
            #mu1=rho*z_now[j]+mu-(-ez[j]+np.sum(z_next*Pi[j,:]))
                Pi[j,:]=int_prob(z_next,mu=mu1,sig=sigma)
            
        
   
       
            
            # if(-mu1+np.sum(z_next*Pi[j,:])<-0.01):
            #     mu1=mu1-(-mu1+np.sum(z_next*Pi[j,:]))
            #     Pi[j,:]=int_prob(z_next,mu=mu1,sig=sigma)
                
            # if(-mu1+np.sum(z_next*Pi[j,:])>0.01):
            #     mu2=mu1+(-mu1+np.sum(z_next*Pi[j,:]))
            #     Pi[j,:]=int_prob(z_next,mu=mu2,sig=sigma)
        

    
    # Pi[:,0] = normcdf_tr( ( z_next[0] + h_half - ez )/sigma)
    # Pi[:,-1] = 1 - normcdf_tr( (z_next[-1] - h_half - ez ) / sigma )
    # for j in range(1,z_next.size - 1):
    #     Pi[:,j] = normcdf_tr( ( z_next[j] + h_half - ez )/sigma) - \
    #         normcdf_tr( ( z_next[j] - h_half - ez )/sigma)
    #for j in range(z_next.size):print(ez[j],np.sum(z_next*Pi[j,:]))
    return Pi
        

def build_s_grid(agrid,n_between,da_min,da_max):
    sgrid = np.array([0.0],agrid.dtype)
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
