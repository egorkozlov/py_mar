#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
aCreated on Tue Sep 17 19:14:08 2019

@author: Egor Kozlov
"""





if __name__ == '__main__':
    
    try:
        from IPython import get_ipython
        get_ipython().magic('reset -f')
    except:
        pass


from platform import system
    
if system() != 'Darwin' and system() != 'Windows':   
    import os
    os.environ['QT_QPA_PLATFORM']='offscreen'



import numpy as np
from residuals import mdl_resid

print('Hi!')




if __name__ == '__main__':
    
    
    

    #Initialize the file with parameters


    x0 = np.array([0.05,0.01,0.02,0.7,0.25])
    
    out, mdl = mdl_resid(x0,return_format=['distance','model'],calibration_report=False,
                         verbose=True)
    print('Done. Residual in point x0 is {}'.format(out))
    
    
   
        

