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
matplotlib.rcParams['axes.unicode_minus'] = False
if __name__ == '__main__':   
    

    
    
    
    def compute(data,alpha=0.05):  
        
        ####################################
        #Regressions
        ####################################
        
        #Do the event study thing around divorce--by wealth
        data=data.dropna()
        ols_mar1 = smf.ols(formula='mar~ C(st)+C(date_rel_month)+C(birth)+female+coll+C(years, Treatment(reference=-3))', 
                          data = data).fit(cov_type='cluster',cov_kwds={'groups': data['st']}) 
        ols_mar2 = smf.ols(formula='mar~ C(st)+C(date_rel_month)+C(birth)+female+coll+C(years, Treatment(reference=-3))', 
                           data = data[data["keep"] ==1]).fit(cov_type='cluster',cov_kwds={'groups': data[data["keep"] ==1]['st']}) 
        ols_mar3 = smf.ols(formula='mar~ C(st)+C(date_rel_month)+C(birth)+female+coll+C(years, Treatment(reference=-3))', 
                           data = data[data["nsfh"] ==1]).fit(cov_type='cluster',cov_kwds={'groups': data[data["nsfh"] ==1]['st']}) 
        ols_mar4 = smf.ols(formula='mar~ C(st)+C(date_rel_month)+C(birth)+female+coll+C(years, Treatment(reference=-3))',
                           data = data[data["nsfh"] ==0]).fit(cov_type='cluster',cov_kwds={'groups': data[data["nsfh"] ==0]['st']}) 
    
        data = data[data["tit"] ==0]
        
        ols_mar1c = smf.ols(formula='mar~ C(st)+C(date_rel_month)+C(birth)+female+coll+C(years, Treatment(reference=-3))', 
                          data = data).fit(cov_type='cluster',cov_kwds={'groups': data['st']}) 
        ols_mar2c = smf.ols(formula='mar~ C(st)+C(date_rel_month)+C(birth)+female+coll+C(years, Treatment(reference=-3))', 
                           data = data[data["keep"] ==1]).fit(cov_type='cluster',cov_kwds={'groups': data[data["keep"] ==1]['st']}) 
        ols_mar3c = smf.ols(formula='mar~ C(st)+C(date_rel_month)+C(birth)+female+coll+C(years, Treatment(reference=-3))', 
                           data = data[data["nsfh"] ==1]).fit(cov_type='cluster',cov_kwds={'groups': data[data["nsfh"] ==1]['st']}) 
        ols_mar4c = smf.ols(formula='mar~ C(st)+C(date_rel_month)+C(birth)+female+coll+C(years, Treatment(reference=-3))',
                           data = data[data["nsfh"] ==0]).fit(cov_type='cluster',cov_kwds={'groups': data[data["nsfh"] ==0]['st']}) 
    
        #Create an Array for the results 
        eventgrid=np.array(np.linspace(-14,10,25),dtype=np.int16) 
        pevent1=np.ones(len(eventgrid))*np.nan 
        pevent2=np.ones(len(eventgrid))*np.nan 
        pevent3=np.ones(len(eventgrid))*np.nan 
        pevent4=np.ones(len(eventgrid))*np.nan 
        pevent1c=np.ones(len(eventgrid))*np.nan 
        pevent2c=np.ones(len(eventgrid))*np.nan 
        pevent3c=np.ones(len(eventgrid))*np.nan 
        pevent4c=np.ones(len(eventgrid))*np.nan 
        pevent1i=np.ones((len(eventgrid),2))*np.nan 
        pevent2i=np.ones((len(eventgrid),2))*np.nan 
        pevent3i=np.ones((len(eventgrid),2))*np.nan 
        pevent4i=np.ones((len(eventgrid),2))*np.nan 
        pevent1ci=np.ones((len(eventgrid),2))*np.nan 
        pevent2ci=np.ones((len(eventgrid),2))*np.nan 
        pevent3ci=np.ones((len(eventgrid),2))*np.nan 
        pevent4ci=np.ones((len(eventgrid),2))*np.nan 
        pevent1i2=np.ones((len(eventgrid),2))*np.nan 
        pevent2i2=np.ones((len(eventgrid),2))*np.nan 
        pevent3i2=np.ones((len(eventgrid),2))*np.nan 
        pevent4i2=np.ones((len(eventgrid),2))*np.nan 
        pevent1ci2=np.ones((len(eventgrid),2))*np.nan 
        pevent2ci2=np.ones((len(eventgrid),2))*np.nan 
        pevent3ci2=np.ones((len(eventgrid),2))*np.nan 
        pevent4ci2=np.ones((len(eventgrid),2))*np.nan 
        
        
        i=0
        for e in eventgrid: 
                   
                if e!=-3.0:
                   
                    try:
                        pevent1[i]=ols_mar1.params['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        pevent2[i]=ols_mar2.params['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]'] 
                        pevent3[i]=ols_mar3.params['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        pevent4[i]=ols_mar4.params['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]'] 
                        pevent1c[i]=ols_mar1c.params['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        pevent2c[i]=ols_mar2c.params['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]'] 
                        pevent3c[i]=ols_mar3c.params['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        pevent4c[i]=ols_mar4c.params['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]'] 
                        pevent1i[i,0]=ols_mar1.conf_int(alpha=alpha)[0]['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        pevent1i[i,1]=ols_mar1.conf_int(alpha=alpha)[1]['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        pevent2i[i,0]=ols_mar2.conf_int(alpha=alpha)[0]['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        pevent2i[i,1]=ols_mar2.conf_int(alpha=alpha)[1]['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        pevent3i[i,0]=ols_mar3.conf_int(alpha=alpha)[0]['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        pevent3i[i,1]=ols_mar3.conf_int(alpha=alpha)[1]['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        pevent4i[i,0]=ols_mar4.conf_int(alpha=alpha)[0]['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        pevent4i[i,1]=ols_mar4.conf_int(alpha=alpha)[1]['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        pevent1ci[i,0]=ols_mar1c.conf_int(alpha=alpha)[0]['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        pevent1ci[i,1]=ols_mar1c.conf_int(alpha=alpha)[1]['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        pevent2ci[i,0]=ols_mar2c.conf_int(alpha=alpha)[0]['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        pevent2ci[i,1]=ols_mar2c.conf_int(alpha=alpha)[1]['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        pevent3ci[i,0]=ols_mar3c.conf_int(alpha=alpha)[0]['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        pevent3ci[i,1]=ols_mar3c.conf_int(alpha=alpha)[1]['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        pevent4ci[i,0]=ols_mar4c.conf_int(alpha=alpha)[0]['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        pevent4ci[i,1]=ols_mar4c.conf_int(alpha=alpha)[1]['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        
                        pevent1i2[i,0]=ols_mar1.conf_int(alpha=0.1)[0]['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        pevent1i2[i,1]=ols_mar1.conf_int(alpha=0.1)[1]['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        pevent2i2[i,0]=ols_mar2.conf_int(alpha=0.1)[0]['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        pevent2i2[i,1]=ols_mar2.conf_int(alpha=0.1)[1]['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        pevent3i2[i,0]=ols_mar3.conf_int(alpha=0.1)[0]['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        pevent3i2[i,1]=ols_mar3.conf_int(alpha=0.1)[1]['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        pevent4i2[i,0]=ols_mar4.conf_int(alpha=0.1)[0]['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        pevent4i2[i,1]=ols_mar4.conf_int(alpha=0.1)[1]['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        pevent1ci2[i,0]=ols_mar1c.conf_int(alpha=0.1)[0]['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        pevent1ci2[i,1]=ols_mar1c.conf_int(alpha=0.1)[1]['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        pevent2ci2[i,0]=ols_mar2c.conf_int(alpha=0.1)[0]['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        pevent2ci2[i,1]=ols_mar2c.conf_int(alpha=0.1)[1]['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        pevent3ci2[i,0]=ols_mar3c.conf_int(alpha=0.1)[0]['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        pevent3ci2[i,1]=ols_mar3c.conf_int(alpha=0.1)[1]['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        pevent4ci2[i,0]=ols_mar4c.conf_int(alpha=0.1)[0]['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']
                        pevent4ci2[i,1]=ols_mar4c.conf_int(alpha=0.1)[1]['C(years, Treatment(reference=-3))[T.'+str(e)+'.0]']

                    except:
                        pevent1[i]=np.nan
                        pevent2[i]=np.nan
                        pevent3[i]=np.nan
                        pevent4[i]=np.nan 
                        pevent1c[i]=np.nan
                        pevent2c[i]=np.nan
                        pevent3c[i]=np.nan
                        pevent4c[i]=np.nan
                        pevent1i[i,0]=np.nan
                        pevent1i[i,1]=np.nan
                        pevent2i[i,0]=np.nan
                        pevent2i[i,1]=np.nan
                        pevent3i[i,0]=np.nan
                        pevent3i[i,1]=np.nan
                        pevent4i[i,0]=np.nan
                        pevent4i[i,1]=np.nan
                        pevent1ci[i,0]=np.nan
                        pevent1ci[i,1]=np.nan
                        pevent2ci[i,0]=np.nan
                        pevent2ci[i,1]=np.nan
                        pevent3ci[i,0]=np.nan
                        pevent3ci[i,1]=np.nan
                        pevent4ci[i,0]=np.nan
                        pevent4ci[i,1]=np.nan
                        
                        pevent1i2[i,0]=np.nan
                        pevent1i2[i,1]=np.nan
                        pevent2i2[i,0]=np.nan
                        pevent2i2[i,1]=np.nan
                        pevent3i2[i,0]=np.nan
                        pevent3i2[i,1]=np.nan
                        pevent4i2[i,0]=np.nan
                        pevent4i2[i,1]=np.nan
                        pevent1ci2[i,0]=np.nan
                        pevent1ci2[i,1]=np.nan
                        pevent2ci2[i,0]=np.nan
                        pevent2ci2[i,1]=np.nan
                        pevent3ci2[i,0]=np.nan
                        pevent3ci2[i,1]=np.nan
                        pevent4ci2[i,0]=np.nan
                        pevent4ci2[i,1]=np.nan
                    i+=1   
                else:
                    pevent1[i]=0.0
                    pevent2[i]=0.0
                    pevent3[i]=0.0
                    pevent4[i]=0.0 
                    pevent1c[i]=0.0
                    pevent2c[i]=0.0
                    pevent3c[i]=0.0
                    pevent4c[i]=0.0
                    pevent1i[i,0]=0.0
                    pevent1i[i,1]=0.0
                    pevent2i[i,0]=0.0
                    pevent2i[i,1]=0.0
                    pevent3i[i,0]=0.0
                    pevent3i[i,1]=0.0
                    pevent4i[i,0]=0.0
                    pevent4i[i,1]=0.0
                    pevent1ci[i,0]=0.0
                    pevent1ci[i,1]=0.0
                    pevent2ci[i,0]=0.0
                    pevent2ci[i,1]=0.0
                    pevent3ci[i,0]=0.0
                    pevent3ci[i,1]=0.0
                    pevent4ci[i,0]=0.0
                    pevent4ci[i,1]=0.0
                    
                    pevent1i2[i,0]=0.0
                    pevent1i2[i,1]=0.0
                    pevent2i2[i,0]=0.0
                    pevent2i2[i,1]=0.0
                    pevent3i2[i,0]=0.0
                    pevent3i2[i,1]=0.0
                    pevent4i2[i,0]=0.0
                    pevent4i2[i,1]=0.0
                    pevent1ci2[i,0]=0.0
                    pevent1ci2[i,1]=0.0
                    pevent2ci2[i,0]=0.0
                    pevent2ci2[i,1]=0.0
                    pevent3ci2[i,0]=0.0
                    pevent3ci2[i,1]=0.0
                    pevent4ci2[i,0]=0.0
                    pevent4ci2[i,1]=0.0
                    i+=1
     
        
       
    
        
        #Return stuff
        return pevent1,pevent2,pevent3,pevent4,pevent1c,pevent2c,pevent3c,pevent4c,pevent1i,pevent2i,pevent3i,pevent4i,pevent1ci2,pevent2ci2,pevent3ci2,pevent4ci2,pevent1i2,pevent2i2,pevent3i2,pevent4i2,pevent1ci2,pevent2ci2,pevent3ci2,pevent4ci2
        
        
        
        
    ############################################
    #Here sample a lot of people from our sample
    ############################################
    
    #vImport the data
    datav=pd.read_stata('NSFG88e.dta')
    datav=datav.drop(['etn','nointc','nch'], axis=1)     
    #datav.dropna(inplace=True)
    
   
    
    #Comupte with true data
    event1,event2,event3,event4,event1c,event2c,event3c,event4c,event1i,event2i,event3i,event4i,event1ci,event2ci,event3ci,event4ci,event1i2,event2i2,event3i2,event4i2,event1ci2,event2ci2,event3ci2,event4ci2=compute(datav)
    
  
    #Create vectors of interst
    eventgrid=np.array(np.linspace(-14,10,25),dtype=np.int16) 
    
    ############################################
    # Event Study 1
    ##########################################   
    fig = plt.figure()    
    f6=fig.add_subplot(1.5,1,1) 
     
    plt.plot(eventgrid, event1,color='r',linestyle='-',  marker='+') 
    plt.fill_between(eventgrid, event1i[:,0], event1i[:,1],alpha=0.2,facecolor='r')
    plt.plot(eventgrid, event1i2[:,0],color='r',linestyle='--',alpha=0.4) 
    plt.plot(eventgrid, event1i2[:,1],color='r',linestyle='--',alpha=0.4) 
    #plt.plot(eventgrid, event1c,color='b',linestyle='-',  marker='x', label='NoTit') 
    #plt.fill_between(eventgrid, event1ci[:,0], event1ci[:,1],alpha=0.2,facecolor='b')
    plt.axvline(x=0.0,color='k')
    plt.plot(eventgrid, event1*0, 'k--') 
    #plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),    
     #         fancybox=True, shadow=True, ncol=2,fontsize=14)    
    plt.xlabel('Event time (years)', fontsize=16)    
    plt.ylabel('Marriage (0/1) wrt baseline', fontsize=16)  
    #plt.savefig('event_1.pgf', bbox_inches = 'tight',pad_inches = 0)  
    #plt.show()
    
    ############################################
    # Event Study 2
    ##########################################   
    fig = plt.figure()    
    f6=fig.add_subplot(1.5,1,1) 
     
    plt.plot(eventgrid, event2,color='r',linestyle='-',  marker='+') 
    plt.fill_between(eventgrid, event2i[:,0], event2i[:,1],alpha=0.2,facecolor='r')
    plt.plot(eventgrid, event2i2[:,0],color='r',linestyle='--',alpha=0.4) 
    plt.plot(eventgrid, event2i2[:,1],color='r',linestyle='--',alpha=0.4) 
    #plt.plot(eventgrid, event2c,color='b',linestyle='-',  marker='x', label='NoTit') 
    #plt.fill_between(eventgrid, event2ci[:,0], event2ci[:,1],alpha=0.2,facecolor='b')
    plt.axvline(x=0.0,color='k')
    plt.plot(eventgrid, event1*0, 'k--') 
    #plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),    
     #         fancybox=True, shadow=True, ncol=2,fontsize=14)    
    plt.xlabel('Event time (years)', fontsize=16)    
    plt.ylabel('Marriage (0/1) wrt baseline', fontsize=16)  
    #plt.savefig('event_2.pgf', bbox_inches = 'tight',pad_inches = 0)  
    #plt.show()
    
    
        ############################################
    # Event Study 3
    ##########################################   
    fig = plt.figure()    
    f6=fig.add_subplot(1.5,1,1) 
     
    plt.plot(eventgrid, event3,color='r',linestyle='-',  marker='+') 
    plt.fill_between(eventgrid, event3i[:,0], event3i[:,1],alpha=0.2,facecolor='r')
    plt.plot(eventgrid, event3i2[:,0],color='r',linestyle='--',alpha=0.4) 
    plt.plot(eventgrid, event3i2[:,1],color='r',linestyle='--',alpha=0.4) 
    #plt.plot(eventgrid, event3c,color='b',linestyle='-',  marker='x', label='NoTit') 
    #plt.fill_between(eventgrid, event3ci[:,0], event3ci[:,1],alpha=0.2,facecolor='b')
    plt.axvline(x=0.0,color='k')
    plt.plot(eventgrid, event1*0, 'k--') 
    #plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),    
     #         fancybox=True, shadow=True, ncol=2,fontsize=14)    
    plt.xlabel('Event time (years)', fontsize=16)    
    plt.ylabel('Marriage (0/1) wrt baseline', fontsize=16)  
    #plt.savefig('event_3.pgf', bbox_inches = 'tight',pad_inches = 0)  
    #plt.show()
    
    
        ############################################
    # Event Study 1
    ##########################################   
    fig = plt.figure()    
    f6=fig.add_subplot(1.5,1,1) 
     
    plt.plot(eventgrid, event4,color='r',linestyle='-',  marker='+') 
    plt.fill_between(eventgrid, event4i[:,0], event4i[:,1],alpha=0.2,facecolor='r')
    plt.plot(eventgrid, event4i2[:,0],color='r',linestyle='--',alpha=0.4) 
    plt.plot(eventgrid, event4i2[:,1],color='r',linestyle='--',alpha=0.4) 
    #plt.plot(eventgrid, event4c,color='b',linestyle='-',  marker='x', label='NoTit') 
    #plt.fill_between(eventgrid, event4ci[:,0], event4ci[:,1],alpha=0.2,facecolor='b')
    plt.axvline(x=0.0,color='k')
    plt.plot(eventgrid, event1*0, 'k--') 
    #plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),    
     #         fancybox=True, shadow=True, ncol=2,fontsize=14)    
    plt.xlabel('Event time (years)', fontsize=16)    
    plt.ylabel('Marriage (0/1) wrt baseline', fontsize=16)  
    #plt.savefig('event_4.pgf', bbox_inches = 'tight',pad_inches = 0)  
    #plt.show()
    
    
    
    #################################
    #Overall Graph of the Effect
    #################################
    fig,axx = plt.subplots()    
    f6=fig.add_subplot(1.5,1,1) 
    
    #Graph for the fit
    datav=datav[datav['nsfh']==1]
    #datav=datav[datav['reln']!=2]
    datav['mar']=datav['mar']*100
    datav['unid']=datav['unid']*100
    ols1=smf.ols(formula='mar~ years',data=datav[(datav['years']>-20) & (datav['years']<=-1)],weights=datav[(datav['years']>-20) & (datav['years']<=-1)]['wgt']).fit()
    ols2=smf.ols(formula='mar~ years',data=datav[(datav['years']>=0) & (datav['years']<=15)],weights=datav[(datav['years']>=0) & (datav['years']<=15)]['wgt']).fit()
    a1=ols1.params['Intercept']
    b1=ols1.params['years']
    a2=ols2.params['Intercept']
    b2=ols2.params['years']
    
    #Graph for the 
    wm = lambda x: np.average(x, weights=data.loc[x.index, "wgt"])
    data=datav[(datav['years']>-20) & (datav['years']<=15)]
    data = data[['mar', 'years']].groupby(datav['years']).agg([wm, 'mean'])
    data.columns = data.columns.droplevel() # remove the multiple levels that were created
    data.columns = ['y_mean', 'y_size', 'x_size', 'x_mean'] # manually set new column names
    
   
    #Scatter
    data.plot.scatter(x='x_mean', y='y_mean', marker='o',color='b',s=30, facecolors='none') # plot
    
    #First
    where=np.array((np.array(data['x_mean']>-20)) & (np.array(data['x_mean']<=-1)),dtype=bool)
    ax=np.array(data['x_mean'])
    plt.plot(ax[where], (ax*b1+a1)[where],color='r',linestyle='--') 
    
    #Second
    where=(np.array(data['x_mean']>=0)) & (np.array(data['x_mean']<=15))
    ax=np.array(data['x_mean'])
    plt.plot(ax[where], (ax*b2+a2)[where],color='r',linestyle='--') 
    
    
    #Horizontal line
    plt.axvline(x=0.0,color='k')
    
    plt.xlabel('Event time (years)', fontsize=18)    
    plt.ylabel(' ') 

   
    #plt.show()
    plt.savefig('event_raw.pgf', bbox_inches = 'tight',pad_inches = 0)  
    #datav.plot.scatter(x='yeard', y='mar')
    
    
       #Graph for the 
    wm = lambda x: np.average(x, weights=data.loc[x.index, "wgt"])
    data=datav[(datav['years']>-20) & (datav['years']<=15)]
    data = data[['mar', 'years']].groupby(datav['years']).agg([wm, 'mean'])
    data.columns = data.columns.droplevel() # remove the multiple levels that were created
    data.columns = ['y_mean', 'y_size', 'x_size', 'x_mean'] # manually set new column names
    
   
    #Scatter
    data.plot.scatter(x='x_mean', y='y_mean', marker='o',color='b',s=30, facecolors='none') # plot
    
    #First
    where=np.array((np.array(data['x_mean']>-20)) & (np.array(data['x_mean']<=-1)),dtype=bool)
    ax=np.array(data['x_mean'])
    plt.plot(ax[where], (ax*b1+a1)[where],color='r',linestyle='--') 
    
    #Second
    where=(np.array(data['x_mean']>=0)) & (np.array(data['x_mean']<=15))
    ax=np.array(data['x_mean'])
    plt.plot(ax[where], (ax*b2+a2)[where],color='r',linestyle='--') 
    
    
    #Horizontal line
    plt.axvline(x=0.0,color='k')
    
    plt.xlabel('event time---years', fontsize=23,labelpad=0)    
    plt.ylabel(' ') 
    #plt.setp(axis='both', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
   
    #plt.show()
    plt.savefig('event_raw_pres.pgf', bbox_inches = 'tight',pad_inches = 0)  
    #datav.plot.scatter(x='yeard', y='mar')
    
    
  