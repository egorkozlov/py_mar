# -*- coding: utf-8 -*-
"""
Decomposition of Welfare Effects

@author: Fabio
"""
import numpy as np
from gridvec import VecOnGrid
import matplotlib.backends.backend_pdf   
from interp_np import interp,interp_np

 
#For nice graphs with matplotlib do the following 
matplotlib.use("pgf") 
matplotlib.rcParams.update({ 
    "pgf.texsystem": "pdflatex", 
    'font.family': 'serif', 
    'font.size' : 11, 
    'text.usetex': True, 
    'pgf.rcfonts': False, 
}) 

def welf_dec(mdl,agents):
    
    ################################################
    #ANALYSIS OF OVERALL WELFARE BY DIVORCE REGIME
    ################################################
 
  

    #Get some stuff from agents
    iexo=agents.iexo     
    state=agents.state   
    theta_t=mdl[0].setup.thetagrid_fine[agents.itheta]   
    setup = mdl[0].setup  
    female=agents.is_female
    cons=agents.c2
    consx=agents.x2
    labor=agents.l2
    labor[state<2]=mdl[0].setup.mlevel
    labori=np.zeros(labor.shape,dtype=np.int32)
    labori[labor>0.5]=1
    psi_check=np.zeros(state.shape)
    betag=mdl[0].setup.pars['beta_t'][0]**(np.linspace(1,len(state[0,:]),len(state[0,:]))-1)
    betam=np.reshape(np.repeat(betag,len(state[:,0])),(len(state[:,0]),len(betag)),order='F')
    
    #Is uni or mut?
    isuni=np.sum(agents.policy_ind,axis=1)==mdl[0].setup.pars['T']
    ismut=np.sum(agents.policy_ind,axis=1)==0
    
    #Fill psi and ushift here
    for i in range(len(state[0,:])):
        psi_check[:,i]=((setup.exogrid.psi_t[i][(setup.all_indices(i,iexo[:,i]))[3]])) 
    
    
    #For welfare-women
    femC=mdl[0].setup.u_part(cons,consx,labori,theta_t,psi_check,mdl[0].setup.pars['u_shift_coh'])[0]*betam
    femM=mdl[0].setup.u_part(cons,consx,labori,theta_t,psi_check,0.0                             )[0]*betam
    femS=mdl[0].setup.u_single_pub(cons,consx,labor)*betam
    
    femU=np.zeros(femC.shape)
    femU[state==3]=femC[state==3]
    femU[state==2]=femM[state==2]
    femU[state==0]=femS[state==0]
    femUt=np.sum(femU,axis=1)
    EF_mut=np.mean(femUt[(female[:,0]==1) & (ismut==1) ])
    EF_uni=np.mean(femUt[(female[:,0]==1) & (isuni==1) ])
    
    #For men
    malC=mdl[0].setup.u_part(cons,consx,labori,theta_t,psi_check,mdl[0].setup.pars['u_shift_coh'])[1]*betam
    malM=mdl[0].setup.u_part(cons,consx,labori,theta_t,psi_check,0.0                             )[1]*betam
    malS=mdl[0].setup.u_single_pub(cons,consx,labor)*betam
    
    malU=np.zeros(malC.shape)
    malU[state==3]=malC[state==3]
    malU[state==2]=malM[state==2]
    malU[state==1]=malS[state==1]
    malUt=np.sum(malU,axis=1)
    EM_mut=np.mean(malUt[(female[:,0]==0) & (ismut==1) ])
    EM_uni=np.mean(malUt[(female[:,0]==0) & (isuni==1) ])#-1.5
    
    ####################################
    #Compute consumption equivalents
    ####################################
    
    #Consumption grid
    consgrid=np.linspace(0.95,1.08,10)
    
    #Util for different cons (women and man)
    valfu=np.ones(consgrid.shape)*100
    valmu=np.ones(consgrid.shape)*100
    valfm=np.ones(consgrid.shape)*100
    valmm=np.ones(consgrid.shape)*100
    
    for i in range(consgrid.size):
        
        #Compute new utility-Women
        femCt=mdl[0].setup.u_part(cons*consgrid[i],consx,labori,theta_t,psi_check,mdl[0].setup.pars['u_shift_coh'])[0]*betam
        femMt=mdl[0].setup.u_part(cons*consgrid[i],consx,labori,theta_t,psi_check,0.0                             )[0]*betam
        femSt=mdl[0].setup.u_single_pub(cons*consgrid[i],consx,labor)*betam
        
        femUt=np.zeros(femCt.shape)
        femUt[state==3]=femCt[state==3]
        femUt[state==2]=femMt[state==2]
        femUt[state==0]=femSt[state==0]
        femUtt=np.sum(femUt,axis=1)
        valfu[i]=np.mean(femUtt[(female[:,0]==1) & (isuni==1) ])
        valfm[i]=np.mean(femUtt[(female[:,0]==1) & (ismut==1) ])
        
        
        #Compute new utility-Men
        malCt=mdl[0].setup.u_part(cons*consgrid[i],consx,labori,theta_t,psi_check,mdl[0].setup.pars['u_shift_coh'])[1]*betam
        malMt=mdl[0].setup.u_part(cons*consgrid[i],consx,labori,theta_t,psi_check,0.0                             )[1]*betam
        malSt=mdl[0].setup.u_single_pub(cons*consgrid[i],consx,labor)*betam
        
        malUt=np.zeros(malCt.shape)
        malUt[state==3]=malCt[state==3]
        malUt[state==2]=malMt[state==2]
        malUt[state==1]=malSt[state==1]
        malUtt=np.sum(malUt,axis=1)
        valmu[i]=np.mean(malUtt[(female[:,0]==0) & (isuni==1) ])#-1.5
        valmm[i]=np.mean(malUtt[(female[:,0]==0) & (ismut==1) ])
        
    #Compute consumption equivalent (using interpolation)
    
    #Women baseline Uni-Mut
    num=np.argmin(~(valfu>EF_mut))
    sh_min=(valfu[num]-EF_mut)/(valfu[num]-valfu[num-1])
    ce_f_u=(sh_min*consgrid[num-1]+(1.0-sh_min)*consgrid[num]-1)*100
   
    #Men baseline Uni-Mut
    num=np.argmin(~(valmu>EM_mut))
    sh_min=(valmu[num]-EM_mut)/(valmu[num]-valmu[num-1])
    ce_m_u=(sh_min*consgrid[num-1]+(1.0-sh_min)*consgrid[num]-1)*100
    
    
    #Computation: CE wrt mutual regime underbaseline
    #to be used for experiments
    
    #Women Uni
    num=np.argmin(~(valfu>-364.40657094625385))
    sh_min=(valfu[num]-(-364.40657094625385))/(valfu[num]-valfu[num-1])
    ce_f_eu=(sh_min*consgrid[num-1]+(1.0-sh_min)*consgrid[num]-1)*100
   
    #Men Uni
    num=np.argmin(~(valmu>-352.28155066914667))
    sh_min=(valmu[num]-(-352.28155066914667))/(valmu[num]-valmu[num-1])
    ce_m_eu=(sh_min*consgrid[num-1]+(1.0-sh_min)*consgrid[num]-1)*100
    
    #Women Mut
    num=np.argmin(~(valfm>-364.40657094625385))
    sh_min=(valfm[num]-(-364.40657094625385))/(valfm[num]-valfm[num-1])
    ce_f_em=(sh_min*consgrid[num-1]+(1.0-sh_min)*consgrid[num]-1)*100
   
    #Men Mut
    num=np.argmin(~(valmm>-352.28155066914667))
    sh_min=(valmm[num]-(-352.28155066914667))/(valmm[num]-valmm[num-1])
    ce_m_em=(sh_min*consgrid[num-1]+(1.0-sh_min)*consgrid[num]-1)*100
    
    print("CE for experiments: f_uni {},m_uni {},f_mut {},m_mut {}"\
          .format(ce_f_eu,ce_m_eu,ce_f_em,ce_m_em))
        
        

    #Create a table with results
    content = r'''\begin{tabular}{cccc}
    \hline\midrule
    \multicolumn{2}{c}{\textbf{Female}}& \multicolumn{2}{c}{\textbf{Male}}\\
    \cmidrule(l){1-2}\cmidrule(l){3-4}
     Mutual Consent & Unilateral Divorce & Mutual Consent & Unilateral Divorce\\
     \cmidrule(l){1-4}
    \multicolumn{4}{c}{\textit{Life-Time utilities in $t=0$}}\\[3ex]
     '''+str(round(EF_mut,2))+''' &'''+str(round(EF_uni,2))+''' &'''+str(round(EM_mut,2))+''' &'''+str(round(EM_uni,2))+''' \\\\
    \cmidrule(l){1-4}
    \multicolumn{4}{c}{\\textit{Welfare Losses with Unilateral Divorce, measured in consumption equivalent variation}}\\\\[3ex]
    \multicolumn{2}{c}{\Chartgirls{'''+str(ce_f_u/10)+'''}}& \multicolumn{2}{c}{\Chartguys{'''+str(ce_m_u/10)+'''}}\\\\[-0.15ex]
    \multicolumn{2}{c}{'''+str(round(ce_f_u,2))+''' \%}& \multicolumn{2}{c}{'''+str(round(ce_m_u,2))+''' \%}\\\\
    \hline\hline
    \end{tabular}
    '''
    
    #Save the File here
    f=open('welfare.tex','w+')
    f.write(content)
    f.close()
    
    # shocks=agents.iexo[:,0] 
    # isfemale=np.array(agents.is_female[:,0],dtype=bool) 
    # ismale=~isfemale 
    # Vf_bil=mdl[0].V[0]['Female, single']['V']
    # Vm_bil=mdl[0].V[0]['Male, single']['V']
    # Vf_uni=mdl[1].V[0]['Female, single']['V'] 
    # Vm_uni=mdl[1].V[0]['Male, single']['V']