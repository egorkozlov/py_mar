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
    

    
    
    
    def compute(data,add1,add2,add3,add4,alpha=0.05):  
        
        ####################################
        #Rich-Poor
        ####################################
        
        #Do the event study thing around divorce--by wealth
        ols_mar = smf.ols(formula='wealth~ C(mary)+C(age)+C(year)+C(less1)+C(event, Treatment(reference=-1))', data = data[data["less1"] ==0]).fit() 
        ols_mar2 = smf.ols(formula='wealth~ C(mary)+C(age)+C(year)+C(less1)+C(event, Treatment(reference=-1))', data = data[data["less1"] ==1]).fit() 
    
    
        #Create an Array for the results 
        eventgrid=np.array(np.linspace(-6,4,6),dtype=np.int16) 
        pevent_mar_wh=np.ones(len(eventgrid))*np.nan 
        pevent_mar_po=np.ones(len(eventgrid))*np.nan 
        pevent_mar_whi=np.ones((len(eventgrid),2))*np.nan 
        pevent_mar_poi=np.ones((len(eventgrid),2))*np.nan 
        
        
        i=0
        for e in eventgrid: 
                   
                if e!=-1.0:
                   
                    try:
                        pevent_mar_wh[i]=ols_mar2.params['C(event, Treatment(reference=-1))[T.'+str(e)+'.0]']
                        pevent_mar_po[i]=ols_mar.params['C(event, Treatment(reference=-1))[T.'+str(e)+'.0]'] 
                        pevent_mar_whi[i,0]=ols_mar2.conf_int(alpha=alpha)[0]['C(event, Treatment(reference=-1))[T.'+str(e)+'.0]']
                        pevent_mar_whi[i,1]=ols_mar2.conf_int(alpha=alpha)[1]['C(event, Treatment(reference=-1))[T.'+str(e)+'.0]']
                        pevent_mar_poi[i,0]= ols_mar.conf_int(alpha=alpha)[0]['C(event, Treatment(reference=-1))[T.'+str(e)+'.0]'] 
                        pevent_mar_poi[i,1]= ols_mar.conf_int(alpha=alpha)[1]['C(event, Treatment(reference=-1))[T.'+str(e)+'.0]'] 
                    except:
                        pevent_mar_wh[i]=np.nan
                        pevent_mar_po[i]=np.nan
                        pevent_mar_whi[i,0]=np.nan
                        pevent_mar_whi[i,1]=np.nan
                        pevent_mar_poi[i,0]=np.nan
                        pevent_mar_poi[i,1]=np.nan
                    
                    i+=1   
                else:
                    pevent_mar_wh[i]=0.0
                    pevent_mar_po[i]=0.0
                    pevent_mar_whi[i,0]=0.0
                    pevent_mar_whi[i,1]=0.0
                    pevent_mar_poi[i,0]=0.0
                    pevent_mar_poi[i,1]=0.0
                    i+=1
                    
        #Sum up rich stuff
        pevent_mar_po=pevent_mar_po+add1
        pevent_mar_wh=pevent_mar_wh+add3
        pevent_mar_poi=pevent_mar_poi+add1
        pevent_mar_whi=pevent_mar_whi+add3
        
       
        ####################################
        #Men-Women
        ####################################
        
        #Do the event study thing around divorce--by wealth
        ols_mar = smf.ols(formula='wealth~C(mary)+ C(age)+C(year)+C(event, Treatment(reference=-1))', data = data[data["sex"] ==1]).fit() 
        ols_mar2 = smf.ols(formula='wealth~C(mary)+ C(age)+C(year)+C(event, Treatment(reference=-1))', data = data[data["sex"] ==2]).fit() 
    
        #Create an Array for the results 
        eventgrid=np.array(np.linspace(-6,4,6),dtype=np.int16) 
        pevent_mar_wo=np.ones(len(eventgrid))*np.nan 
        pevent_mar_me=np.ones(len(eventgrid))*np.nan 
        pevent_mar_woi=np.ones((len(eventgrid),2))*np.nan 
        pevent_mar_mei=np.ones((len(eventgrid),2))*np.nan 
        
        
        i=0
        for e in eventgrid: 
                   
                if e!=-1.0:
                   
                    try:
                        pevent_mar_wo[i]=ols_mar.params['C(event, Treatment(reference=-1))[T.'+str(e)+'.0]']
                        pevent_mar_me[i]=ols_mar2.params['C(event, Treatment(reference=-1))[T.'+str(e)+'.0]'] 
                        pevent_mar_woi[i,0]=ols_mar.conf_int(alpha=alpha)[0]['C(event, Treatment(reference=-1))[T.'+str(e)+'.0]']
                        pevent_mar_woi[i,1]=ols_mar.conf_int(alpha=alpha)[1]['C(event, Treatment(reference=-1))[T.'+str(e)+'.0]']
                        pevent_mar_mei[i,0]= ols_mar2.conf_int(alpha=alpha)[0]['C(event, Treatment(reference=-1))[T.'+str(e)+'.0]'] 
                        pevent_mar_mei[i,1]= ols_mar2.conf_int(alpha=alpha)[1]['C(event, Treatment(reference=-1))[T.'+str(e)+'.0]'] 
                    
                    except:
                        
                        pevent_mar_wh[i]=np.nan
                        pevent_mar_po[i]=np.nan
                        pevent_mar_woi[i,0]=np.nan
                        pevent_mar_woi[i,1]=np.nan
                        pevent_mar_mei[i,0]=np.nan
                        pevent_mar_mei[i,1]=np.nan
                    
                    i+=1   
                else:
                    pevent_mar_wo[i]=0.0
                    pevent_mar_me[i]=0.0
                    pevent_mar_woi[i,0]=0.0
                    pevent_mar_woi[i,1]=0.0
                    pevent_mar_mei[i,0]=0.0
                    pevent_mar_mei[i,1]=0.0
                    i+=1
                    
        #Sum up rich stuff
        pevent_mar_me=pevent_mar_me+add4
        pevent_mar_wo=pevent_mar_wo+add2
        pevent_mar_mei=pevent_mar_mei+add4
        pevent_mar_woi=pevent_mar_woi+add2
        
        
        #Return stuff
        return pevent_mar_wh,pevent_mar_po,pevent_mar_me,pevent_mar_wo,pevent_mar_whi,pevent_mar_poi,pevent_mar_mei,pevent_mar_woi
        
        
        
        
    ############################################
    #Here sample a lot of people from our sample
    ############################################
    
    #vImport the data
    datav=pd.read_stata('sepe.dta')
    #datav.dropna(inplace=True)
    
    addo1=np.mean(datav.loc[(datav['event']==-2) & (datav['less1']==0),'wealth'])
    addo2=np.mean(datav.loc[(datav['event']==-2) & (datav['sex']==1),'wealth'])
    addo3=np.mean(datav.loc[(datav['event']==-2) & (datav['less1']==1),'wealth'])
    addo4=np.mean(datav.loc[(datav['event']==-2) & (datav['sex']==2),'wealth'])
    
    
    #Comupte with true data
    event_mar_wh,event_mar_po,event_mar_me,event_mar_wo,event_mar_whi,event_mar_poi,event_mar_mei,event_mar_woi=compute(datav,addo1,addo2,addo3,addo4,alpha=0.05)
    
  
    #Create vectors of interst
    eventgrid=np.array(np.linspace(-6,4,6),dtype=np.int16) 
    
    ############################################
    # Event Study Love Shock 
    ##########################################   
    fig = plt.figure()    
    f6=fig.add_subplot(1.5,1,1) 
     
    plt.plot(eventgrid, event_mar_po/1000,color='r',linestyle='-',  marker='+', label='Poor') 
    plt.fill_between(eventgrid, event_mar_poi[:,0]/1000, event_mar_poi[:,1]/1000,alpha=0.2,facecolor='r')
    plt.plot(eventgrid, event_mar_wh/1000,color='b',linestyle='--', marker='x', label='Rich') 
    plt.fill_between(eventgrid, event_mar_whi[:,0]/1000, event_mar_whi[:,1]/1000,alpha=0.2,facecolor='b')
    plt.axvline(x=0.0,color='k')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),    
              fancybox=True, shadow=True, ncol=2,fontsize=14)    
    plt.xlabel('Event time (Years)', fontsize=16)    
    plt.ylabel('Net Worth (\$ 1000s)', fontsize=16)  
    plt.savefig('eventwh_r_p_c.pgf', bbox_inches = 'tight',pad_inches = 0)  
    #plt.show()
    
    ##########################################    
    # Event Study Love Shock 
    ##########################################   
    fig = plt.figure()    
    f6=fig.add_subplot(1.5,1,1) 
     
    plt.plot(eventgrid, event_mar_wo/1000,color='r',linestyle='-',  marker='+', label='Women') 
    plt.fill_between(eventgrid, event_mar_woi[:,0]/1000, event_mar_woi[:,1]/1000,alpha=0.2,facecolor='r')
    plt.plot(eventgrid, event_mar_me/1000,color='b',linestyle='--', marker='x', label='Men') 
    plt.fill_between(eventgrid, event_mar_mei[:,0]/1000, event_mar_mei[:,1]/1000,alpha=0.2,facecolor='b')
    plt.axvline(x=0.0,color='k')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),    
              fancybox=True, shadow=True, ncol=2, fontsize=14)    
    plt.xlabel('Event time (Years)', fontsize=16)    
    plt.ylabel('Net Worth (\$ 1000s)', fontsize=16)  
    plt.savefig('eventwh_r_s_c.pgf', bbox_inches = 'tight',pad_inches = 0)  
    #plt.show()
    
    