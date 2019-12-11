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

def find_between( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""


#f = line.split()


from p_function import fun


try:
    mkdir('Job')
except:    
    pass
    
    
start = default_timer()
tic = default_timer()


time_abort = 1200.0 # time when the worker stops, in seconds



while True:
    toc = default_timer()
    if toc - start > time_abort: break
        
    
    if toc - tic > 5.0:
        print('I am ok, running for {:.1f} min'.format((toc-start)/60))
        tic = toc
        
    
    sleep(rs()) # sleep random number of seconds
    
    
    
    li_txt = [f for f in listdir('Job') if f.endswith('.txt') and f.startswith('in')]
        
    if len(li_txt) == 0: continue

    fname = li_txt[0]
    
    try:
        num = int(find_between(fname,'in','.txt'))
    except:
        num = 0
        print('something wrong with file named {}'.format(fname))
    
    
    fname_full = 'Job/{}'.format(fname)
    
    file_in = open(fname_full)
    try:        
        st = file_in.readline().split()
        x = [float(s) for s in st]
        file_in.close()
        remove(fname_full)
        print('I got a job to solve {}'.format(fname))
        f = fun(x)
    except KeyboardInterrupt:
        raise KeyboardInterrupt()
    except BaseException as e:
        print('error while solving {}'.format(fname))
        print('error text is ' + str(e))
        #try:
        #    remove(fname)
        #except:
        #    pass
        continue
        f = 1e6          
    
    
    try:
        file_tmp_name = 'Job/tmp{}.txt'.format(num)
        file_tmp = open(file_tmp_name,'w+')  

        try: # if f is list
            x_str = [str(x) for x in f]
            f = ' '.join(x_str)
        except:
            pass
        
        file_tmp.write('{}'.format(f))
        file_tmp.close()        
        file_out_name = 'Job/out{}.txt'.format(num)
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
        
    
chdir(cwd)
