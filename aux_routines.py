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



import cupy as cp
def cp_take_along_axis(a, indices, axis):
    """Take values from the input array by matching 1d index and data slices.
    Args:
        a (cupy.ndarray): Array to extract elements.
        indices (cupy.ndarray): Indices to take along each 1d slice of ``a``.
        axis (int): The axis to take 1d slices along.
    Returns:
        cupy.ndarray: The indexed result.
    .. seealso:: :func:`numpy.take_along_axis`
    """

    if indices.dtype.kind not in ('i', 'u'):
        raise IndexError('`indices` must be an integer array')

    if axis is None:
        a = a.ravel()
        axis = 0

    ndim = a.ndim

    if not (-ndim <= axis < ndim):
        raise _errors._AxisError('Axis overrun')

    axis %= a.ndim

    if ndim != indices.ndim:
        raise ValueError(
            '`indices` and `a` must have the same number of dimensions')

    fancy_index = []
    for i, n in enumerate(a.shape):
        if i == axis:
            fancy_index.append(indices)
        else:
            ind_shape = (1,) * i + (-1,) + (1,) * (ndim - i - 1)
            fancy_index.append(cp.arange(n).reshape(ind_shape))

    return a[fancy_index]

