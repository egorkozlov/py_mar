# -*- coding: utf-8 -*-
"""
Implementation of TikTak code as described in:
    
    'Benchmarking Global Optimizers'
     by Antoine Arnoud,Fatih Guvenen and Tatjana Kleineberg'

@author: Fabio
"""
import sobol_seq
import numpy as np
from scipy.optimize import minimize

def tiktak(nthreads,N,N_st,xl,xu,f,tole):
    
    #Initial cheks
    assert len(xl)==len(xu)
    
    assert N>=N_st
    
    #First Create a Sobol Sequence
    init = sobol_seq.i4_sobol_generate(len(xl),N) # generate many draws from uniform
    #init=init[:,0]   
    
    #Get point on the grid
    x_init=xl*(1-init)+xu*init
    x_init=x_init.T

    #Get fitness of initial points
    fx_init=np.ones(N)
    for j in range(N):
        print(j,555555555555)
        fx_init[j]=f(x_init[:,j])
    
    #Sort in ascending order of fitness
    order=np.argsort(fx_init,axis=0)
    fx_init=fx_init[order]
    x_init=x_init[:,order]
    
    #Take only the first N_st realization
    fx_init=fx_init[0:N_st]
    x_init=x_init[:,0:N_st]
   
    
    #List containing parameters
    param=list()
    
    #Initialize the Iteration
    for i in range(N_st):
        print(99999999999,i)
        #Sort lists
        def sortFirst(val): 
           return val[0]
        
        #Get the starting point for local minimization
        if i>0:
            
            #Open file with all paramters:
            #file = open('wisdom.dat','r')
            param.sort(key=sortFirst)
            print('f best so far is {} and x is {}'.format(param[0][0],param[0][1]))
            xm=param[0][1]
            dump=min(max(0.1,((i+1)/N_st)**(0.5)),0.995)
            xc=dump*xm+(1-dump)*x_init[:,i]
        else:
            xc=x_init[:,0]
            
        #Local Minimization
        res = minimize(f,xc,method='Nelder-Mead',tol=tole,options={'maxiter':10})
        print('Final value is {}'.format(f(res.x)))
        
        #Store the new value in wisdom
        #file = open('wisdom.dat',"w")
        #file.write(np.stack)
        #file.close()
        param=param+[(res.fun,res.x)]
        
    param.sort(key=sortFirst)
    return param[0]
    
#if __name__ == '__main__':
    
 #   def ff(x):
  #     return x[0]**2/200+x[1]**2/200+x[2]**2/200-np.cos(x[2]/3**(0.5))*np.cos(x[0])*np.cos(x[1]/2**(0.5))+2#(x[0]-1)**2+(x[1]-2)**2
    
  #  param=tiktak(1,30,21,np.array([-10,-14,-30]),np.array([21,35,3]),ff,1e-4)