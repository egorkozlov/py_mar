#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 16:50:59 2020

@author: egorkozlov
"""
import numpy as np
from scipy.optimize import brentq

def int_sol(m_in,newton=True,step=1e-6,nint=2000,*,A,alp,sig,xi,lam,kap,lbr):
    m_in = np.atleast_1d(m_in)    
    
    def foc_expression(x):
        return ((A/alp)**(1/sig))*(x**lam + kap*(1-lbr)**lam)**((lam+xi-1)/(lam*sig)) * \
                x**((1-lam)/sig) + x
    
    
    
    rf_res = lambda x : foc_expression(x) - m_in[0]
    
    x0 = 0
    
    for _ in range(5):
        try:
            x0 = brentq(rf_res,step,m_in[0]-step)
            break
        except:
            step = step/10
            pass
        
        
    
    xgrid = np.linspace(x0,m_in.max(),nint)
    
    
    m_implied = foc_expression(xgrid)
          
    x_interpolated = np.maximum(np.interp(m_in,m_implied,xgrid),x0)
    
    if not newton:
        xout = x_interpolated
    else:
        # perform 1 iteration of Newton method
        def foc_deriv(x):
            logder = ((lam+xi-1)/(sig)) * (x**(lam-1))/(x**lam + kap*(1-lbr)**lam) + \
                ((1-lam)/sig)/x
            return 1 + (foc_expression(x) - x)*logder
    
        f_res = foc_expression(x_interpolated) - m_in
        f_der = foc_deriv(x_interpolated)
    
        x_improved = np.maximum( x_interpolated - f_res/f_der, x0)
    
        x_interpolated = x_improved
            
        xout = x_interpolated
    
    def util(x,M):
        c = M - x
        uc = A*c**(1-sig)/(1-sig)
        assert np.all(uc < 0)
        ux = alp*(x**lam + kap*(1-lbr)**lam)**((1-xi)/lam)/(1-xi)
        assert np.all(ux < 0)
        return uc + ux
    
    cout = m_in - xout
    uout = util(xout,m_in)
    
    return xout, cout, uout