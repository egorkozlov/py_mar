#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This collects auxiliary functions

"""
import numpy as np



def num_to_nparray(*things_in):
    # makes everything numpy array even if passed as number
    
    out = ( 
            x if isinstance(x,np.ndarray) else np.array([x]) 
            for x in things_in
          ) # this is generator not tuple
    return tuple(out)


def unify_sizes(*things_in):
    # makes all the inputs the same size. basically broadcasts all singletones
    # to be the same size as the first non-singleton array. Also checks is all
    # non-singleton arrays are of the same size
    
    shp = 1
    for x in things_in:
        if x.size != 1: 
            shp = x.shape
            break
        
    if shp != 1:
        assert all([(x.size==1 or x.shape==shp) for x in things_in]), 'Wrong shape, it is {}, but we broadcast it to {}'.format([x.shape for x in things_in],shp)
        things_out = tuple( ( x if x.shape==shp else x*np.ones(shp) for x in things_in))
    else:
        things_out = things_in
        
    return things_out




def first_true(mask, axis=None, invalid_val=-1):
    return np.where(mask.any(axis=axis), mask.argmax(axis=axis), invalid_val)

def last_true(mask, axis=None, invalid_val=-1):
    val = mask.shape[axis] - np.flip(mask, axis=axis).argmax(axis=axis) - 1
    return np.where(mask.any(axis=axis), val, invalid_val)


