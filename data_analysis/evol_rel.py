# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 20:00:31 2020

Regression+nice graphs for paper on cohabitation and uni divorce

@author: Fabio
"""
import pandas as pd
import numpy as np   
import pickle    
import statsmodels.formula.api as smf 
import matplotlib.pyplot as plt  
import matplotlib.backends.backend_pdf   


 #For nice graphs with matplotlib do the following  
matplotlib.use("pgf")  
matplotlib.rcParams.update({  
    "pgf.texsystem": "pdflatex",  
    'font.family': 'sans-serif',  
    'font.size' : 14,  
    'text.usetex': True,  
    'pgf.rcfonts': False,  
})  

if __name__ == '__main__':   

    
    datav=pd.read_stata('NSFG88e2.dta')
    datav['mar']=datav['mar']*100
    datav['unid']=datav['unid']*100
    datav=datav[datav['nsfh']==1]
    #datav=datav[datav['reln']!=2]
    
    # Define a lambda function to compute the weighted mean:
    wm = lambda x: np.average(x, weights=datav.loc[x.index, "wgt"])
    
    #################################
    #Graphs for evolution over time
    #################################
    fig = plt.figure()    
    f6=fig.add_subplot(1.5,1,1) 
    
     #Graph for the 
     #data=datav[(datav['date_rel']>=-20) & (datav['years']<=15)]
    data=datav[(datav['date_rel']>=1955) & (datav['date_rel']<=1988)]
    
    data = data[['mar', 'date_rel','unid']].groupby(datav['date_rel']).agg([wm, 'size'])
    data.columns = data.columns.droplevel() # remove the multiple levels that were created
    data.columns = ['y_mean', 'y_size', 'x_mean', 'x_size','z_mean','z_size'] # manually set new column names
    ax=np.array(data['x_mean'])
    yx=np.array(data['y_mean'])
    zx=np.array(data['z_mean'])
    
    fig, ax2 = plt.subplots()

    color = 'b'
    ax2.set_xlabel('Year the relationship started', fontsize=18)
    #ax2.set_ylabel('% new couples choosing marriage', color=color, fontsize=18)  # we already handled the x-label with ax1
    plt.plot(ax,yx, marker='x',color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    #plt.setp(ax2.get_xticklabels(), fontsize=14)
    #plt.setp(ax2.get_yticklabels(), fontsize=14)
    fig.tight_layout()  # otherwise the right y-label is slightly clipped

    
    ax1 = ax2.twinx()  # instazzzzntiate a second axes that shares the same x-axis
    
    color = 'r'
    
    #ax1.set_ylabel('% relationships born under U.D.', color=color, fontsize=18)
    plt.plot(ax,zx, marker='+', color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    #plt.setp(ax1.get_xticklabels(), fontsize=14)
    #plt.setp(ax1.get_yticklabels(), fontsize=14)
   
    
    
    #plt.xlabel('Year', fontsize=16)    
    #plt.ylabel('% couples that marry', fontsize=16) 
   
   
    plt.savefig('years_raw.pgf', bbox_inches = 'tight',pad_inches = 0)  
    
    
    #For presentation: before...
    fig, ax2 = plt.subplots()

    color = 'b'
    ax2.set_xlabel('year the relationship started', fontsize=22,labelpad=0)
    #ax2.set_ylabel('% new couples choosing marriage', color=color, fontsize=18)  # we already handled the x-label with ax1
    plt.plot(ax,yx, marker='x',color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    plt.setp(ax2.get_xticklabels(), fontsize=20)
    plt.setp(ax2.get_yticklabels(), fontsize=20)
    fig.tight_layout()  # otherwise the right y-label is slightly clipped

    
    ax1 = ax2.twinx()  # instazzzzntiate a second axes that shares the same x-axis
    
    color = 'white'
    
    plt.plot(ax,zx, marker='+', color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    plt.setp(ax1.get_xticklabels(), fontsize=20)
    plt.setp(ax1.get_yticklabels(), fontsize=20)
   
    
    
    #plt.xlabel('Year', fontsize=16)    
    #plt.ylabel('% couples that marry', fontsize=16) 
   
   
    plt.savefig('years_raw_before.pgf', bbox_inches = 'tight',pad_inches = 0)  
    
    fig, ax2 = plt.subplots()

    color = 'b'
    ax2.set_xlabel('year the relationship started', fontsize=22,labelpad=0)
    #ax2.set_ylabel('% new couples choosing marriage', color=color, fontsize=18)  # we already handled the x-label with ax1
    plt.plot(ax,yx, marker='x',color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    plt.setp(ax2.get_xticklabels(), fontsize=20)
    plt.setp(ax2.get_yticklabels(), fontsize=20)
    fig.tight_layout()  # otherwise the right y-label is slightly clipped

    
    ax1 = ax2.twinx()  # instazzzzntiate a second axes that shares the same x-axis
    
    color = 'red'
    
    plt.plot(ax,zx, marker='+', color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    plt.setp(ax1.get_xticklabels(), fontsize=20)
    plt.setp(ax1.get_yticklabels(), fontsize=20)
   
    
    
    #plt.xlabel('Year', fontsize=16)    
    #plt.ylabel('% couples that marry', fontsize=16) 
   
   
    plt.savefig('years_raw_after.pgf', bbox_inches = 'tight',pad_inches = 0)  
    
    
    #Graph of difference by education with data from Manning
    fig = plt.figure()    
    f6=fig.add_subplot(1.5,1,1) 
     
    time=np.array([1987,1995,2002,2010])
    edu=np.array([31,37,45,50])
    nedu=np.array([32,50,63,67])
    plt.plot(time,edu,color='b',  marker='+') 
    plt.plot(time, nedu,color='r',marker='x') 
    #plt.axvline(x=0.0,color='k')
    #plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),    
     #         fancybox=True, shadow=True, ncol=2, fontsize=14)    
    #plt.xlabel('Interview Years', fontsize=16)    
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlabel('interview year', fontsize=22)  
    plt.savefig('cbye.pgf', bbox_inches = 'tight',pad_inches = 0) 