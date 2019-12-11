#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 17:07:47 2019

@author: egorkozlov
"""

import os
from time import sleep
from timeit import default_timer



def fun_check(x):
    #sleep(5.0)
    return sum([i**2 for i in x]), x[0]

#f = line.split()



def compute_for_values(values,timeout=240.0,print_every=15.0):
    cwd = os.getcwd()
    
    if cwd.endswith('Job'): os.chdir('..')   
 
    
    try:
        os.chdir('Job')            
    except:
        os.mkdir('Job')
        os.chdir('Job')
        
    
    
    assert type(values) is list
    
    names_in = list()
    names_out = list()
    
    # find names of in and out
    # clear files if something already exists
    for ival in range(len(values)):
        namein = 'in{}.txt'.format(ival)
        
        try:
            os.remove(namein)
        except:
            pass
        
        names_in.append(namein)
        nameout = 'out{}.txt'.format(ival)
        names_out.append(nameout)
        
        try:
            os.remove(nameout)
        except:
            pass
        
    
        
        
    def create_in(fname,x):
        x_str = [str(y) for y in x]
        towrite = ' '.join(x_str)
        file_in = open(fname,'w+')  
        file_in.write(towrite)
        file_in.close()
    
    
    # create bunch of in files and write values of them
    for fname, x in zip(names_in,values):
        create_in(fname,x)
        
        
        
    
    time_start = [0.0]*len(names_in)
    time_took = [0.0]*len(names_in)
    started = [False]*len(names_in)
    finished = [False]*len(names_in)
    
    start = default_timer()
    tic = default_timer()
    
    
    while True:
        
        sleep(1.0)
        
        
        ld = os.listdir()
        
        
        li_in =  [f for f in ld if f.endswith('.txt') and f.startswith('in') ]
        li_out = [f for f in ld if f.endswith('.txt') and f.startswith('out')]
         
        for i, name in enumerate(names_in):
            if (name not in li_in): 
                if not started[i]:
                    started[i] = True # mark as started
                    time_start[i] = default_timer()
                elif not finished[i]:
                    # check if takes too long
                    time_since = default_timer() - time_start[i]
                    
                    if time_since > timeout: # if does restart
                        print('timeout for i = {}, recreating'.format(i))
                        time_start[i] = 0
                        started[i] = False
                        create_in(name,values[i])
                        
            if (names_out[i] in li_out) and (not finished[i]) and (started[i]):
                finished[i] = True
                time_took[i] = default_timer() - time_start[i]
                

        if set(li_out) == set(names_out):
            break # job is done!
            
        # time stats  sometimes if not done
        toc = default_timer()
        
        if toc - tic > print_every:
            print('{} started, {} finished, {} not started, running for {:.1f} minutes'.
                  format(sum(started),sum(finished),len(values)-sum(started),(toc-start)/60))
            tic = toc
            
            
    
    
    fout = list()
    for i, name in enumerate(names_out):
        file = open(name)
        
        # this handles both lists and values
        val = [float(x) for x in file.readline().split()]
        if len(val)==1: val = val[0]            
        
        fout.append(val)
        file.close()
        os.remove(name)
        
        
        
    os.chdir(cwd)
    return fout
    
if __name__ == '__main__':
    vals = [[1,2],[4.0,5.0,6.0],[3.0],[-1,-2,-3]]
    vout = compute_for_values(vals)
    print(vout)
 
    
    