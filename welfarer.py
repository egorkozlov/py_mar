# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 10:41:24 2020

@author: Fabio
"""

import numpy as np   
import matplotlib.pyplot as plt   
#from matplotlib.pyplot import plot, draw, show   
import matplotlib.backends.backend_pdf   
from statutils import strata_sample
import os,glob,subprocess
 
#For nice graphs with matplotlib do the following 
matplotlib.use("pgf") 
matplotlib.rcParams.update({ 
    "pgf.texsystem": "pdflatex", 
    'font.family': 'serif', 
    'font.size' : 11, 
    'text.usetex': True, 
    'pgf.rcfonts': False, 
}) 

def welfare(mdl,agents):
    
    #Welafare of Agents in t=0 under the two different regimes  
    shocks=agents.iexo[:,0]
    isfemale=np.array(agents.is_female[:,0],dtype=bool)
    ismale=~isfemale
    Vf_bil=mdl[0].V[0]['Female, single']['V'][0,shocks[isfemale]]
    Vm_bil=mdl[0].V[0]['Male, single']['V'][0,shocks[ismale]]
    Vf_uni=mdl[1].V[0]['Female, single']['V'][0,shocks[isfemale]]
    Vm_uni=mdl[1].V[0]['Male, single']['V'][0,shocks[ismale]]
    changef=(Vf_bil-Vf_uni)/np.abs(Vf_bil)
    changem=(Vm_bil-Vm_uni)/np.abs(Vm_bil) 
    
    
    #Write a latex table that contains those Values
    content = r'''\begin{tabular}{@{} l l l @{}}
	Pull Up Method & \Chartguys{1.000} & \Chartguys{0.600}\\
	Move Field     & \Chartguys{0.269} & \Chartguys{0.783}
    \end{tabular}
    '''
    
    content = r'''\begin{tabular}{lcccc}
    \hline\hline
    &\multicolumn{2}{c}{\textbf{Female}}& \multicolumn{2}{c}{\textbf{Male}}\\ \cmidrule(l){2-3}\cmidrule(l){4-5}
    & \textit{Mutual Consent} & \textit{Unilateral Divorce} & \textit{Mutual Consent} & \textit{Unilateral Divorce}\\
	& \Chartgirls{'''+str(np.mean(Vf_bil))+'''} & \Chartgirls{'''+str(np.mean(Vf_uni))+'''} & \Chartguys{'''+str(np.mean(Vm_bil))+'''} & \Chartguys{'''+str(np.mean(Vm_uni))+'''}\\
    Utility & '''+str(np.mean(Vf_bil))+''' &'''+str(np.mean(Vf_uni))+''' &'''+str(np.mean(Vm_bil))+''' &'''+str(np.mean(Vm_uni))+''' \\ 
    \% change & \multicolumn{2}{c}{'''+str(np.mean(changem))+'''} & \multicolumn{2}{c}{'''+str(np.mean(changef))+'''} \\
    \hline\hline
    \end{tabular}
    '''
 
    #Save the File here
    f=open('welfare.tex','w+')
    f.write(content)
    f.close()
    
#    commandLine = subprocess.Popen(['pdflatex', 'welfare.tex'])
#    commandLine.communicate()
#
#    os.unlink('welfare.tex')
#    os.unlink('welfare.log')
#    os.unlink('welfare.aux')