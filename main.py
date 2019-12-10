#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
aCreated on Tue Sep 17 19:14:08 2019

@author: Egor Kozlov
"""


if __name__ == '__main__':
    
    
    #Clean Memory
    
    try:
        from IPython import get_ipython
        get_ipython().magic('reset -f')
    except:
        pass
    
    #If on server set Display
    from platform import system
    
    if system() != 'Darwin' and system() != 'Windows':   
        import os
        os.environ['QT_QPA_PLATFORM']='offscreen'
 
    import numpy as np
    from model import Model
    from setup import DivorceCosts
  
    

            
    #Create grids of parameters
    sigma_psi_g=np.linspace(0.01,0.3,3)
    sigma_psi_init_g=np.linspace(0.05,0.5,3)
    di_co_g=np.linspace(0.05,0.3,3)
    bila=np.array([False,True])
    
    
    
    #Initialize the file with parameters

    import xlwt 
    from xlwt import Workbook 
   
    #File that gives some info on the go
    f = open("iterations","a").close()
    # Workbook is created 
    wb = Workbook() 
    
    # add_sheet is used to create sheet. 
    sheet1 = wb.add_sheet('Sheet 1', cell_overwrite_ok=True) 
    sheet1.write(0, 0, 'Bilateral Divorce')
    sheet1.write(0, 1, 'sigma_psi') 
    sheet1.write(0, 2, 'sigma_psi_init') 
    sheet1.write(0, 3, 'u_cost') 
    sheet1.write(0, 4, '% coh before ret')
    sheet1.write(0, 5, '% mar berfore ret')
    sheet1.write(0, 6, 'hazd[average]')
    sheet1.write(0, 7, 'hazm[average]')
    sheet1.write(0, 8, 'hazs[average]')
    sheet1.write(0, 9, 'flsm')
    sheet1.write(0, 10, 'flsc')
    sheet1.write(0, 11, '% coh mean')
    sheet1.write(0, 12, '% mar mean')
    
    
    
    
    row=0
    for i in range(len(sigma_psi_g)):
        for j in range(len(sigma_psi_init_g)):
            for k in range(len(di_co_g)):
                for h in bila:
                
                

                    row=row+1
                    
                    f = open("iterations.txt","w")
                    f.write('{}'.format(row))
                    dc = DivorceCosts(unilateral_divorce=h,assets_kept = 1.0,u_lost_m=di_co_g[k],u_lost_f=di_co_g[k],eq_split=0.0)
                    sc = DivorceCosts(unilateral_divorce=True,assets_kept = 1.0,u_lost_m=0.00,u_lost_f=0.00)
                    mdl = Model(iterator_name='default',
                                divorce_costs=dc,separation_costs=sc,sigma_psi=sigma_psi_g[i],sigma_psi_init=sigma_psi_init_g[j])
                    
                    graphs=True
                    #gassets,iexo,state,gtheta=mdl.solve_sim()
                    mdl.solve_sim(simulate=True)
                    #gassets, iexo, state, gtheta = mdl.agents.gsavings_c, mdl.agents.iexo, mdl.agents.state, mdl.agents.gtheta
                    
                    Tret = mdl.setup.pars['Tret']
                    
                    #Write results on spreadsheet
                    sheet1.write(row, 0, '{}'.format(h)) 
                    sheet1.write(row, 1, '{}'.format(sigma_psi_g[i])) 
                    sheet1.write(row, 2, '{}'.format(sigma_psi_init_g[j])) 
                    sheet1.write(row, 3, '{}'.format(di_co_g[k])) 
                    sheet1.write(row, 4, '{}'.format(mdl.moments['share coh'][Tret-1])) 
                    sheet1.write(row, 5, '{}'.format(mdl.moments['share mar'][Tret-1])) 
                    sheet1.write(row, 6, '{}'.format(np.mean(mdl.moments['hazard div']))) 
                    sheet1.write(row, 7, '{}'.format(np.mean(mdl.moments['hazard mar']))) 
                    sheet1.write(row, 8, '{}'.format(np.mean(mdl.moments['hazard sep'])))
                    sheet1.write(row, 9, '{}'.format(np.mean(mdl.moments['flsm'][1:-1])))
                    sheet1.write(row, 10, '{}'.format(np.mean(mdl.moments['flsc'][1:-1])))
                    sheet1.write(row, 11, '{}'.format(np.mean(mdl.moments['share coh'][:Tret]))) 
                    sheet1.write(row, 12, '{}'.format(np.mean(mdl.moments['share mar'][:Tret]))) 
                

               
    wb.save("model_parameters.xls") 
    f.close()           
    #Graphs Here
    
    
    #Indexes for the graphs
    if graphs:
        ai=0
        zfi=0
        zmi=4
        psii=5
        ti=1
        thi=10
        
        #Actual Graphs
        mdl.graph(ai,zfi,zmi,psii,ti,thi)
        
        #If you plan to use graphs only once, deselect below to save space on disk
        #os.remove('name_model.pkl') 
    


   
        


    