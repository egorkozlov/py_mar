#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 17:07:47 2019

@author: egorkozlov
"""
#from time import sleep
from main import mdl_resid


# we need fun() to be possible, do not remove None
def fun(x=None):
    if x is None:
        return mdl_resid()
    else:
        return mdl_resid(x)
    
    