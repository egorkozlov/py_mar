# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 16:32:47 2020

@author: Fabio

Statistical Utilities for Structural Estimation

"""

def strata_sample(namevar,dsample,frac=0.1,weights=False):
    
    #Genrate array for pop count
    namevar_1=[None]*len(namevar)
    number=[None]*len(namevar)
    num=0
    for i in namevar:
        
        namevar_1[num]=namevar[num][1:-1]
        number[num]=str(num)
        num+=1
        
    #Genrate distribution
    pop_count = dsample.groupby(namevar_1)[namevar_1[0]].count()
   
    #Preparation for actual sampling
    subset=''
    for r,i in zip(namevar,number):
        subset=subset+'(     population[{}] == pop_count.index[x][{}] )&'.format(r,i) 
        
    #Get weights right
    if not weights:
        dsample['weights']=1.0
        
    #Stratify Sample
    
    def map_fun(x):
        p = pop_count
        return dsample[eval(subset[0:-1])].sample(frac=frac,weights=dsample['weights'])
    
    stratified_sample= list(map(map_fun, range(len(pop_count))))

    
    return pd.concat(stratified_sample) 
    




if __name__ == '__main__':   
       
    import pandas as pd
    import numpy as np
    
    # Generate random population (100K)
    
    population = pd.DataFrame(index=range(0,100000))
    population['income'] = 0
    population['income'].iloc[39000:80000] = 1
    population['income'].iloc[80000:] = 2
    population['sex'] = np.random.randint(0,2,100000)
    population['age'] = np.random.randint(0,4,100000)
    
    #pop_count = population.groupby(['income', 'sex', 'age'])['income'].count()
    
    

    
    sample2=strata_sample(["'income'", "'sex'", "'age'"],population,frac=0.1)
    