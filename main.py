#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
aCreated on Tue Sep 17 19:14:08 2019

@author: Egor Kozlov
"""

#import numpy as np
from model import Model
import os



if __name__ == '__main__':
    
    #Solve the model
    mdl = Model(iterator_name='default-timed',
                divorce_costs={'unilateral_divorce':True})
    
    graphs=True
    #gassets,iexo,state,gtheta=mdl.solve_sim()
    mdl.solve_sim(simulate=True)
    #gassets, iexo, state, gtheta = mdl.agents.gsavings_c, mdl.agents.iexo, mdl.agents.state, mdl.agents.gtheta
    mdl.time_statistics()
    
    
    #Graphs Here
    
    
    #Indexes for the graphs
    if graphs:
        ai=0
        zfi=0
        zmi=4
        psii=5
        ti=0
        thi=10
        
        #Actual Graphs
        Packed,dec=mdl.graph(ai,zfi,zmi,psii,ti,thi)
        
        #If you plan to use graphs only once, deselect below to save space on disk
        #os.remove('name_model.pkl') 
    
   
    
    