#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 09:13:02 2020

@author: egorkozlov
"""

from model import Model
import numpy as np


x0 = np.array([0.2,0.0710307,1.11501,0.543047,0.050264,0.005,-0.09])

from calibration_params import calibration_params

pars = calibration_params()[-1](x0)
print(pars)

pars.pop('alost')

mdl = Model(**pars,verbose=True,solve_till=-3)


decisions = mdl.decisions[-2]['Female, single']['Decision']
