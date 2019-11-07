#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 20:13:57 2019

@author: egorkozlov
"""

from interp_np import interp            

# this is griddend linear interpolant that uses interpolation routines from
# interp_np. It can be applied to arrays specifying values of function to get
# their interpolated values along an arbitrary dimesnion (see method apply)

class VecOnGrid(object):
    def __init__(self,grid,values,trim=True):
        self.val = values
        self.grid = grid
        self.i, self.wnext = interp(self.grid,self.val,return_wnext=True,trim=trim)
        self.n = self.i.size
        self.wthis = 1-self.wnext
        self.trim = trim
        
    def apply(self,xin,axis=0,take=None,pick=None,reshape_i=True):
        # this applies interpolator to array xin along dimension axis
        # and additionally takes indices specified in take list
        # take's elements are assumed to be 2-element tuple where take[0] 
        # is axis and take[1] is indices. 
        
        
        nd = xin.ndim
        assert axis < nd
        
        if isinstance(take,tuple):
            take = [take]
        
        ithis = [slice(None)]*nd
        inext = [slice(None)]*nd
        
        
        
        
        if pick is None:        
            ithis[axis] = self.i
            inext[axis] = self.i + 1
            wthis = self.wthis
            wnext = self.wnext
            n = self.n
        else:
            ithis[axis] = self.i[pick]
            inext[axis] = self.i[pick] + 1
            wthis = self.wthis[pick]
            wnext = self.wnext[pick]
            n = pick.size
            
            
        
        shp_i = (1,)*axis + (n,) + (1-axis)*(1,)
        # TODO: this is not general but let's see if we need more
        # this creates 2-dimensional thing
        shp_w = (1,)*axis + (n,) + (nd-1-axis)*(1,)
        
        
        if reshape_i:
            ithis[axis] = ithis[axis].reshape(shp_i)
            inext[axis] = inext[axis].reshape(shp_i)
            
        
        if take is not None:
            
            dimextra = [t[0] for t in take]
            iextra = [t[1] for t in take]
            for d, i in zip(dimextra,iextra):
                ithis[d] = i
                inext[d] = i
                
                
            if not reshape_i: # remove one dimension from w
                # if reshape_i is false it assumes that each index in interpolant
                # corresponds to each index in take[1], and therefore the
                # dimensionality is reduced accordingly (I know this seems hard to 
                # comprehend but something like V[[0,1],[0,1],:] returns 2-dimensional
                # array, and V[ [[0],[1]], [[0,1]], :] returns 3-dimensional array, so
                # reshape_i=False leads to dimensionality reduction 
                #assert all([i.shape == self.i.shape for i in iextra])
                shp_w = list(shp_w)
                shp_w = [s for i, s in enumerate(shp_w) if i not in dimextra] # remove extra elements
                shp_w = tuple(shp_w)
                
                
        wthis = wthis.reshape(shp_w)
        wnext = wnext.reshape(shp_w)
        
        
        return wthis*xin[tuple(ithis)] + wnext*xin[tuple(inext)]
    
    
    def update(self,where,values):
        # this replaces values at positions where with values passed to the
        # function and recomputes interpolation weights
        
        assert where.size == values.size
        
        self.val[where] = values
        self.i[where], self.wnext[where] = \
            interp(self.grid,self.val[where],return_wnext=True,trim=self.trim)
        self.wthis[where] = 1 - self.wnext[where]
        