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
import dfogn
import dfols
from time import sleep
import pickle

def tiktak(nthreads,N,N_st,xl,xu,f,tole=1e-3,nelder=True,refine=False):
    
    #Initial cheks
    assert len(xl)==len(xu)
    
    assert N>=N_st
    
    #Redesign the function
    q= lambda x:np.array([f(x)])
    
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
    fx_init=np.ones(N)
    for j in range(N):
       
        fx_init[j]=f(x_init[:,j])
    
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
    param=list()
    filer('wisdom.pkl',param,True)
         
    #Initialize loop value
    i=-1
    filer('loop.pkl',i,True)
    
    ############################################
    #2 LOCAL STAGE
    ###########################################
    #ite=0
    #Check if stop
    while i<=(N_st-2):
       
        #Update Iteration
        i=filer('loop.pkl',i,False)
        
        i+=1
        
        #Write on file
        filer('loop.pkl',i,True)

        
        #Sort lists
        def sortFirst(val): 
           return val[0]
        
        #Get the starting point for local minimization
        if i>0:
            
            #Open File with best solution so far
            param=filer('wisdom.pkl',param,False)
                 
            param.sort(key=sortFirst)
            print('f best so far is {} and x is {}'.format(param[0][0],param[0][1]))
            xm=param[0][1]
            
            #Get right sobol sequence point
            xt=filer('sobol.pkl',i,False)
            
            #Determine the initial position
            dump=min(max(0.1,((i+1)/N_st)**(0.5)),0.995)
            
            xc=dump*xm+(1-dump)*xt[:,i]
            
        else:
            xc=x_init[:,0]
            
        #Local Minimization
        
        if nelder:
            
            res = minimize(q,xc,method='Nelder-Mead',tol=tole)#,options={'maxiter':300}
            #ite=ite+res.nfev        
            print('Final value is {}'.format(f(res.x)))
            param=param+[(res.fun,res.x)]
        else:

            #res=dfogn.solve(q, xc,rhoend=tole)
            res=dfols.solve(q, xc,rhoend=tole)
            #ite=ite+res.nf
            print('Final value is {}'.format(f(res.x)))
            param=param+[(res.f,res.x)]

        
        #Save Updated File
        filer('wisdom.pkl',param,True)
        
        param.sort(key=sortFirst)
        
            
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
                with open(filename, 'rb+') as file:
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