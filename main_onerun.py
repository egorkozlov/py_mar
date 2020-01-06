
1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
18
19
20
21
22
23
24
25
26
27
28
29
30
31
32
33
34
35
36
37
38
39
40
41
42
43
44
45
46
47
48
49
50
51
52
53
54
55
56
57
58
59
60
61
62
63
64
65
66
67
68
69
70
71
72
73
74
75
76
77
78
	
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
from data_moments import dat_moments
import pickle
 
print('Hi!')
 
 
 
 
if __name__ == '__main__':
     
     
     
    #Build  data moments and pickle them
    packed_stuff=dat_moments(100,weighting=False)
    with open('moments.pkl', 'wb+') as file:
        pickle.dump(packed_stuff,file)
         
    #Initialize the file with parameters
 
 
    x0 = np.array([0.071875 ,0.006875,0.02125,0.59375,0.171875])#p.exp(np.array([ -1.8603,-8.1430,-1.57934,0.25130,-0.4991]))#0.08512367 -0.03874894 -0.05721577  0.57536013  0.20720013
    #0.08512367,-0.03874894,-0.05721577, 0.57536013, 0.20720013
 
 
    out, mdl = mdl_resid(x0,return_format=['distance','model'],calibration_report=False,
                         verbose=True)
    print('Done. Residual in point x0 is {}'.format(out))
     
    graphs=True
    #Indexes for the graphs
    if graphs:
        ai=30
        zfi=0
        zmi=4
        psii=3
        ti=1
        thi=10
         
        #Actual Graphs
        #mdl.graph(ai,zfi,zmi,psii,ti,thi)
         
        #If you plan to use graphs only once, deselect below to save space on disk
        #os.remove('name_model.pkl')
     
     
    
        
