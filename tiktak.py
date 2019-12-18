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
from p_client import compute_for_values
from time import sleep
import pickle

def tiktak(nthreads,N,N_st,xl,xu,f,tole=1e-3,nelder=True,refine=False):
    
    #Initial cheks
    assert len(xl)==len(xu)
    
    assert N>=N_st
    
    
    ############################################
    #1 INITIALIZATION
    ###########################################
    
    #First Create a Sobol Sequence
    init = sobol_seq.i4_sobol_generate(len(xl),N) # generate many draws from uniform
    #init=init[:,0]   
    
    #Get point on the grid
    x_init=xl*(1-init)+xu*init
    x_init=x_init.T
    x_init=x_init.squeeze()

    #Get fitness of initial points
    
    pts = [ ('compute',x_init[:,j]) for j in range(N)]
    fx_init = compute_for_values(pts)
    fx_init = np.array(fx_init)
    #Sort in ascending order of fitness
    order=np.argsort(fx_init,axis=0)
    fx_init=fx_init[order]
    x_init=x_init[:,order]
    
    #Take only the first N_st realization
    fx_init=fx_init[0:N_st]
    x_init=x_init[:,0:N_st]
   
    #Create a file with sobol sequence points
    filer('sobol.pkl',x_init,True)    
    
    #List containing parameters and save them in file
    param=list([ (fx_init[0], x_init[:,0])])
    filer('wisdom.pkl',param,True)
         
    
    vals = [('minimize',(i,N_st)) for i in range(N_st)]
    
    compute_for_values(vals)
    
    param = filer('wisdom.pkl',None,write=False)
    
    ############################################
    #3 TOPPING RULE
    ###########################################
    #print(999,ite)
    #Final Refinement
    if refine:
        res = minimize(f,param[0][1],method='Nelder-Mead',tol=1e-8)
        param[0]=(res.fun,res.x)
    
    
    return param[0]
    
##########################################
#Functions
#########################################
    
#Write on Function
def filer(filename,array,write=True):
    
    while True:
        try:
            if write:
                with open(filename, 'wb+') as file:
                    pickle.dump(array,file)
            else:
                with open(filename, 'rb') as file:
                    array=pickle.load(file)
                return array
                
            break
        except KeyboardInterrupt:
            raise KeyboardInterrupt()
        except:
            print('Problems opening the file {}'.format(filename))
            sleep(0.5)
    

##########################################
# UNCOMMENT BELOW FOR TESTING
######################################### 
#def ff(x):
#    return 10*3+1+(x[0]**2-10*np.cos(2*np.pi*x[0]))+(x[1]**2-10*np.cos(2*np.pi*x[1]))+(x[2]**2-10*np.cos(2*np.pi*x[2]))
#param=tiktak(1,100,30,np.array([-25.12,-7.12,-5.12]),np.array([15.12,50.12,1.12]),ff,1e-3,nelder=False,refine=False)