#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 17:07:47 2019

@author: egorkozlov
"""
#from time import sleep
from main import mdl_resid



# we need fun() to be possible, do not remove None
def fun(x):
    assert type(x) is tuple, 'x must be a tuple!'
    
    action = x[0]
    args = x[1]
    
    assert type(action) is str, 'x[0] should be string for action'
    assert len(x) <= 2, 'too many things in x! x is (action,agrs)'
    
    
    if action == 'test':
        return mdl_resid()
    elif action == 'compute':
        return mdl_resid(args)
    elif action == 'minimize':	
        from scipy.optimize import minimize
        x0 = args
        print('I am running a local maximizer starting at {}'.format(x0))
        res = minimize(mdl_resid,x0,method='Nelder-Mead',options={'maxiter':20})
        return res.x, res.fun
    else:
        raise Exception('unsupported action or format')
    
    
