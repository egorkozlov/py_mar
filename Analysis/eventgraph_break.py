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
    'font.family': 'serif',  
    'font.size' : 11,  
    'text.usetex': True,  
    'pgf.rcfonts': False,  
})  

if __name__ == '__main__':   
    

    
    
    
    def compute(data,add1,add2):  
        
        ####################################
        #Rich-Poor
        ####################################
        
        #Do the event study thing around divorce--by wealth
        ols_mar = smf.ols(formula='wealth2~ C(dif1)+C(age)+C(year)+C(less1)+C(less1)*C(event, Treatment(reference=-1))', data = data).fit() 
    
        #Create an Array for the results 
        eventgrid=np.array(np.linspace(-6,6,7),dtype=np.int16) 
        pevent_mar_wh=np.ones(len(eventgrid))*np.nan 
        pevent_mar_po=np.ones(len(eventgrid))*np.nan 
        
        
        i=0
        for e in eventgrid: 
                   
                if e!=-2.0:
                   
                    try:
                        pevent_mar_wh[i]=ols_mar.params['C(less1)[T.1.0]:C(event, Treatment(reference=-1))[T.'+str(e)+'.0]']
                        pevent_mar_po[i]=ols_mar.params['C(event, Treatment(reference=-1))[T.'+str(e)+'.0]'] 
                    except:
                        pevent_mar_wh[i]=np.nan
                        pevent_mar_po[i]=np.nan
                    
                    i+=1
                else:
                    pevent_mar_wh[i]=0.0
                    pevent_mar_po[i]=0.0
                    i+=1
                    
        #Sum up rich stuff
        pevent_mar_po=pevent_mar_po+add1
        pevent_mar_wh=pevent_mar_wh+pevent_mar_po+ols_mar.params['C(less1)[T.1.0]']
        
       
        ####################################
        #Men-Women
        ####################################
        
        #Do the event study thing around divorce--by wealth
        ols_mar = smf.ols(formula='wealth2~C(dif1)+ C(age)+C(year)+C(sex, Treatment(reference=1))+C(sex, Treatment(reference=1))*C(event, Treatment(reference=-1))', data = data).fit() 
    
        #Create an Array for the results 
        eventgrid=np.array(np.linspace(-6,6,7),dtype=np.int16) 
        pevent_mar_wo=np.ones(len(eventgrid))*np.nan 
        pevent_mar_me=np.ones(len(eventgrid))*np.nan 
        
        
        i=0
        for e in eventgrid: 
                   
                if e!=-2.0:
                   
                    try:
                        pevent_mar_wo[i]=ols_mar.params['C(sex, Treatment(reference=1))[T.2.0]:C(event, Treatment(reference=-1))[T.'+str(e)+'.0]']
                        pevent_mar_me[i]=ols_mar.params['C(event, Treatment(reference=-1))[T.'+str(e)+'.0]'] 
                    
                    except:
                        pevent_mar_wh[i]=np.nan
                        pevent_mar_po[i]=np.nan
                    
                    i+=1
                else:
                    pevent_mar_wo[i]=0.0
                    pevent_mar_me[i]=0.0
                    i+=1
                    
        #Sum up rich stuff
        pevent_mar_me=pevent_mar_me+add2
        pevent_mar_wo=pevent_mar_wo+pevent_mar_me+ols_mar.params['C(sex, Treatment(reference=1))[T.2.0]']

        
        
        #Return stuff
        return pevent_mar_wh,pevent_mar_po,pevent_mar_me,pevent_mar_wo
        
        
        
        
    ############################################
    #Here sample a lot of people from our sample
    ############################################
    
    #vImport the data
    datav=pd.read_stata('sepe.dta')
    #datav.dropna(inplace=True)
    
    addo1=np.mean(datav.loc[(datav['event']==-2) & (datav['less1']==0),'wealth'])
    addo2=np.mean(datav.loc[(datav['event']==-2) & (datav['sex']==1),'wealth'])
    
    #Comupte with true data
    event_mar_wh,event_mar_po,event_mar_me,event_mar_wo=compute(datav,addo1,addo2)
    
    #Set number of repetitions 
    boot=100

    
    
    #We draw peopple, so I subset the data such data only on observation is left
    datad=datav.loc[datav['id1']==1,'id']
    n=len(datad)    
    nn=n*boot 
    data_big=datad.sample(n=nn,replace=True,random_state=4) 
    
    
    #Create vectors of interst
    eventgrid=np.array(np.linspace(-6,6,7),dtype=np.int16) 
    event_mar_poB=np.zeros((len(eventgrid),boot))   
    event_mar_whB=np.zeros((len(eventgrid),boot))   
    event_mar_woB=np.zeros((len(eventgrid),boot))   
    event_mar_meB=np.zeros((len(eventgrid),boot))  
    

    
    for i in range(boot):   
       
        #First take the needed indexes
        a1t=data_big[(i*n):((i+1)*n)].copy().reset_index()   
        
        #Then connect with the full dataset
        a1=pd.merge(datav, a1t, on='id')
        
        #Finally compute stuff
        event_mar_whB[:,i],event_mar_poB[:,i],event_mar_meB[:,i],event_mar_woB[:,i]=compute(a1.copy(),addo1,addo2)   

    #Get the confidence intervals
    event_mar_whi=np.array((np.percentile(event_mar_whB,2.5,axis=1),np.percentile(event_mar_whB,97.5,axis=1)))   
    event_mar_poi=np.array((np.percentile(event_mar_poB,2.5,axis=1),np.percentile(event_mar_poB,97.5,axis=1)))   
    event_mar_mei=np.array((np.percentile(event_mar_meB,2.5,axis=1),np.percentile(event_mar_meB,97.5,axis=1)))   
    event_mar_woi=np.array((np.percentile(event_mar_woB,2.5,axis=1),np.percentile(event_mar_woB,97.5,axis=1)))   
    
    
    
    ############################################
    # Event Study Love Shock 
    ##########################################   
    fig = plt.figure()    
    f6=fig.add_subplot(1.5,1,1) 
     
    plt.plot(eventgrid, event_mar_po,color='r',linestyle='-',  marker='+', label='Poor') 
    plt.fill_between(eventgrid, event_mar_poi[0,:], event_mar_poi[1,:],alpha=0.2,facecolor='r')
    plt.plot(eventgrid, event_mar_wh,color='b',linestyle='--', marker='x', label='Rich') 
    plt.fill_between(eventgrid, event_mar_whi[0,:], event_mar_whi[1,:],alpha=0.2,facecolor='b')
    plt.axvline(x=0.0,color='k')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),    
              fancybox=True, shadow=True, ncol=2,fontsize=14)    
    plt.xlabel('Event time (Years)', fontsize=16)    
    plt.ylabel('Net Worth', fontsize=16)  
    plt.savefig('eventwh_r_p_c.pgf', bbox_inches = 'tight',pad_inches = 0)  
    #plt.show()
    
    ##########################################    
    # Event Study Love Shock 
    ##########################################   
    fig = plt.figure()    
    f6=fig.add_subplot(1.5,1,1) 
     
    plt.plot(eventgrid, event_mar_wo,color='r',linestyle='-',  marker='+', label='Women') 
    plt.fill_between(eventgrid, event_mar_woi[0,:], event_mar_woi[1,:],alpha=0.2,facecolor='r')
    plt.plot(eventgrid, event_mar_me,color='b',linestyle='--', marker='x', label='Men') 
    plt.fill_between(eventgrid, event_mar_mei[0,:], event_mar_mei[1,:],alpha=0.2,facecolor='b')
    plt.axvline(x=0.0,color='k')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),    
              fancybox=True, shadow=True, ncol=2, fontsize=14)    
    plt.xlabel('Event time (Years)', fontsize=16)    
    plt.ylabel('Net Worth', fontsize=16)  
    plt.savefig('eventwh_r_s_c.pgf', bbox_inches = 'tight',pad_inches = 0)  
    #plt.show()
    
    