#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 17:07:47 2019

@author: egorkozlov
"""
#from time import sleep
from main import mdl_resid

def fun(x):
    #return (x[0]-1)**2        
    return mdl_resid(x)