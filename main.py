#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 19:14:08 2019

@author: Egor Kozlov
"""

#import numpy as np
from model import Model



            


if __name__ == '__main__':
    
    mdl = Model(iterator_name='default-timed',
                divorce_costs={'unilateral_divorce':True})
    mdl.solve_sim()
    mdl.time_statistics()
    