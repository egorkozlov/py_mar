#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 17:07:47 2019

@author: egorkozlov
"""
import numpy as np
from time import sleep

def fun(x):
    sleep(5.0)
    return np.array([sum([i**2 for i in x]), x[0]])