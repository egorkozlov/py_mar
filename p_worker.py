#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 17:07:47 2019

@author: egorkozlov
"""

from os import listdir, remove, getcwd, chdir, mkdir
from shutil import copyfile
from time import sleep
from timeit import default_timer
from numpy.random import random_sample as rs

import pickle

def find_between( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""


#f = line.split()


from p_function import fun

# this tests if the worker can run at least something
try:
    out = fun()
    print('I was able to get the value {} in the default point'.format(out))
except:
    raise Exception('I cannot compute the function in the default point, sorry...')


try:
    mkdir('Job')
except:    
    pass
    
    
start = default_timer()
tic = default_timer()


time_abort = 7200.0 # time when the worker stops, in seconds



while True:
    toc = default_timer()
    if toc - start > time_abort: break
        
    
    if toc - tic > 5.0:
        print('I am ok, running for {:.1f} min'.format((toc-start)/60))
        tic = toc
        
    
    sleep(rs()) # sleep random number of seconds
    
    
    
    li_txt = [f for f in listdir('Job') if f.endswith('.pkl') and f.startswith('in')]
        
    if len(li_txt) == 0: continue

    fname = li_txt[0]
    
    try:
        num = int(find_between(fname,'in','.pkl'))
    except:
        num = 0
        print('something wrong with file named {}'.format(fname))
    
    
    fname_full = 'Job/{}'.format(fname)
    
    file_in = open(fname_full,'rb')
    try:                
        x = pickle.load(file_in)
        file_in.close()
        remove(fname_full)
        print('I got a job to solve {}'.format(fname))
        f = fun(x)
    except KeyboardInterrupt:
        raise KeyboardInterrupt()
    except BaseException as e:
        print('error while solving {}'.format(fname))
        print('error text is ' + str(e))        
        continue
        f = 1e6          
    
    
    try:
        file_tmp_name = 'Job/tmp{}.pkl'.format(num)
        file_tmp = open(file_tmp_name,'wb+')  
        pickle.dump(f,file_tmp)
        file_tmp.close()        
        file_out_name = 'Job/out{}.pkl'.format(num)
        try:
            remove(file_out_name)
        except:
            pass
        
        copyfile(file_tmp_name,file_out_name)
        remove(file_tmp_name)   
    except:
        try:
            file_tmp.close()
        except:
            pass
        
        try:
            remove(file_tmp_name)
        except:
            pass
        
        print('could not get output!')
        continue
        
