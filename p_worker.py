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


debug_mode = True

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
    print('I am trying to test myself...')
    out = fun(('test',0))
    print('I was able to get the value {} in the default point'.format(out))
except KeyboardInterrupt:
    raise KeyboardInterrupt
except BaseException as a:
    print('I cannot compute the function in the default point, sorry...')
    raise a
   


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
        try:
            x = pickle.load(file_in)
        except EOFError:
            print('file {} seems empty... deleting it just in case'.format(fname_full))
            file_in.close()
            remove(fname_full)
            continue
        
        file_in.close()
        remove(fname_full)
        print('I got a job to solve {}'.format(fname))
        print('request is {}'.format(x))
        
        try:
            f = fun(x)
        except KeyboardInterrupt:
            raise KeyboardInterrupt()
        except BaseException as e:
            print('error while evaluating function at {}'.format(fname))
            print('error text is ' + str(e)) 
            if debug_mode: raise e
            f = 1e6
            
    except KeyboardInterrupt:
        raise KeyboardInterrupt()
    except BaseException as e:
        print('error while solving {}'.format(fname))
        print('error text is ' + str(e))
        if debug_mode: raise e        
        continue
    
    
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
        
