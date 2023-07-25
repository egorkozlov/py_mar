"""
Graphs for paper on cohabitation and uni divorce for Figure 1
@author: Fabio
"""

import pandas as pd
import numpy as np     
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
matplotlib.rcParams['axes.unicode_minus'] = False


if __name__ == '__main__':   
    
    #Path where to save output
    mypath=r'C:\Users\Fabio\Dropbox\py_mar_3\Analysis\output'
    
    ##########################################    
    ##########################################
    #####   PANEL A) OF FIGURE 1
    #########################################
    ##########################################
    datav=pd.read_stata(mypath+'/fig1.dta')
    datav['mar']=datav['mar']*100
    datav['unid']=datav['unidd']*100
    wm = lambda x: np.average(x, weights=datav.loc[x.index, "wgt"])
    datav['date_rel']=datav['date_rel_agg']+0.0
 
    #Prepre and subset the data 
    data=datav[(datav['date_rel']>=1950)]
    data = data[['mar', 'date_rel','unid']].groupby(datav['date_rel']).agg([wm, 'size'])
    data.columns = data.columns.droplevel() # remove the multiple levels that were created
    data.columns = ['y_mean', 'y_size', 'x_mean', 'x_size','z_mean','z_size'] # manually set new column names
    ax=np.array(data['x_mean'])
    yx=np.array(data['y_mean'])
    zx=np.array(data['z_mean'])
    
      
    ########################################################
    #Figure 1 for the paper (also to use in presentation)
    ########################################################
    fig, ax2 = plt.subplots()    
    color = 'b'
    ax2.set_xlabel('Year the relationship started', fontsize=18)
    plt.plot(ax,yx, marker='x',color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    fig.tight_layout() 
    ax1 = ax2.twinx() 
    color = 'r'
    plt.plot(ax,zx, marker='+', color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    
    #Save
    plt.savefig(mypath+'/years_raw.pgf', bbox_inches = 'tight',pad_inches = 0) 
    plt.savefig(mypath+'/years_raw.png', bbox_inches = 'tight',pad_inches = 0) 
    
    
    ########################################################
    #Figure 1 for the presentation only (only one line)
    ########################################################
    fig, ax2 = plt.subplots()  
    color = 'b'
    ax2.set_xlabel('Year the relationship started', fontsize=18)
    plt.plot(ax,yx, marker='x',color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    fig.tight_layout()  
    ax1 = ax2.twinx()    
    color = 'white'   
    plt.plot(ax,zx, marker='+', color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    
    #Save
    plt.savefig(mypath+'/years_raw_before.pgf', bbox_inches = 'tight',pad_inches = 0)  
    plt.savefig(mypath+'/years_raw_before.png', bbox_inches = 'tight',pad_inches = 0)  
    
 
      
    ##########################################    
    ##########################################
    #####   PANEL B) OF FIGURE 1
    #########################################
    ##########################################

    
    #################################
    #Overall Graph of the Effect
    #################################
    fig,axx = plt.subplots()    
    f6=fig.add_subplot(1.5,1,1) 
    
    #Prepare and subset the data
    data=datav[(datav['years']>-20) & (datav['years']<=15)]
    data = data[['mar', 'years']].groupby(datav['years']).agg([wm, 'mean'])
    data.columns = data.columns.droplevel() # remove the multiple levels that were created
    data.columns = ['y_mean', 'y_size', 'x_size', 'x_mean'] # manually set new column names
    
    
    #Regressions to be added on the graph
    ols1=smf.ols(formula='mar~ years',data=datav[(datav['years']>-20) & (datav['years']<=-1)],weights=datav[(datav['years']>-20) & (datav['years']<=-1)]['wgt']).fit()
    ols2=smf.ols(formula='mar~ years',data=datav[(datav['years']>=0) & (datav['years']<=15)],weights=datav[(datav['years']>=0) & (datav['years']<=15)]['wgt']).fit()
    a1=ols1.params['Intercept']
    b1=ols1.params['years']
    a2=ols2.params['Intercept']
    b2=ols2.params['years']
    
   
    ###############################
    #Graphs elements below
    ##############################
    
    #Scatter
    data.plot.scatter(x='x_mean', y='y_mean', marker='o',color='b',s=30, facecolors='none') # plot
    
    #First regression
    where=np.array((np.array(data['x_mean']>-20)) & (np.array(data['x_mean']<=-1)),dtype=bool)
    ax=np.array(data['x_mean'])
    plt.plot(ax[where], (ax*b1+a1)[where],color='r',linestyle='--') 
    
    #Second regression
    where=(np.array(data['x_mean']>=0)) & (np.array(data['x_mean']<=15))
    ax=np.array(data['x_mean'])
    plt.plot(ax[where], (ax*b2+a2)[where],color='r',linestyle='--') 
    
    #Horizontal line
    plt.axvline(x=0.0,color='k')
    
    #Labels
    plt.xlabel('Event time (years)', fontsize=18)    
    plt.ylabel(' ') 
    

    #Save
    plt.savefig(mypath+'/event_raw.pgf', bbox_inches = 'tight',pad_inches = 0) 
    plt.savefig(mypath+'/event_raw.png', bbox_inches = 'tight',pad_inches = 0) 
    
   