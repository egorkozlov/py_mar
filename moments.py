
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
79
80
81
82
83
84
85
86
87
88
89
90
91
92
93
94
95
96
97
98
99
100
101
102
103
104
105
106
107
108
109
110
111
112
113
114
115
116
117
118
119
120
121
122
123
124
125
126
127
128
129
130
131
132
133
134
135
136
137
138
139
140
141
142
143
144
145
146
147
148
149
150
151
152
153
154
155
156
157
158
159
160
161
162
163
164
165
166
167
168
169
170
171
172
173
174
175
176
177
178
179
180
181
182
183
184
185
186
187
188
189
190
191
192
193
194
195
196
197
198
199
200
201
202
203
204
205
206
207
208
209
210
211
212
213
214
215
216
217
218
219
220
221
222
223
224
225
226
227
228
229
230
231
232
233
234
235
236
237
238
239
240
241
242
243
244
245
246
247
248
249
250
251
252
253
254
255
256
257
258
259
260
261
262
263
264
265
266
267
268
269
270
271
272
273
274
275
276
277
278
279
280
281
282
283
284
285
286
287
288
289
290
291
292
293
294
295
296
297
298
299
300
301
302
303
304
305
306
307
308
309
310
311
312
313
314
315
316
317
318
319
320
321
322
323
324
325
326
327
328
329
330
331
332
333
334
335
336
337
338
339
340
341
342
343
344
345
346
347
348
349
350
351
352
353
354
355
356
357
358
359
360
361
362
363
364
365
366
367
368
369
370
371
372
373
374
375
376
377
378
379
380
381
382
383
384
385
386
387
388
389
390
391
392
393
394
395
396
397
398
399
400
401
402
403
404
405
406
407
408
409
	
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 10:26:49 2019
 
This file comupte simulated moments + optionally
plots some graphs
 
@author: Fabio
"""
 
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.pyplot import plot, draw, show
import matplotlib.backends.backend_pdf
import pickle
 
def moment(mdl,draw=True):
#This function compute moments coming from the simulation
#Optionally it can also plot graphs about them. It is feeded with
#matrixes coming from simulations
 
    agents = mdl.agents
    #Import simulated values
    assets_t=agents.setup.agrid_c[agents.iassets]
    iexo=agents.iexo
    state=agents.state
    theta_t=agents.setup.thetagrid_fine[agents.itheta]
    setup = agents.setup
     
    mdl.moments = dict()
     
     
    #As a first thing we unpack assets and theta
    N=len(state)
     
    #Get states codes
    state_codes = {name: i for i, name in enumerate(agents.setup.state_names)}
     
    ###########################################
    #Moments: FLS over time by Relationship
    ###########################################
     
     
    flsm=np.ones(agents.setup.pars['T'])
    flsc=np.ones(agents.setup.pars['T'])
     
     
    for t in range(agents.setup.pars['T']):
         
        pick = agents.state[:,t]==2       
        if pick.any(): flsm[t] = np.array(setup.ls_levels)[agents.ils_i[pick,t]].mean()
        pick = agents.state[:,t]==3
        if pick.any(): flsc[t] = np.array(setup.ls_levels)[agents.ils_i[pick,t]].mean()
         
     
         
    mdl.moments['flsm'] = flsm
    mdl.moments['flsc'] = flsc
     
    ###########################################
    #Moments: Variables over Age
    ###########################################
     
    relt=np.zeros((len(state_codes),agents.setup.pars['T']))
    ass_rel=np.zeros((len(state_codes),agents.setup.pars['T']))
    inc_rel=np.zeros((len(state_codes),agents.setup.pars['T']))
     
     
     
    for ist,sname in enumerate(state_codes):
        for t in range(agents.setup.pars['T']):
             
             
            ftrend = agents.setup.pars['f_wage_trend'][t]
            mtrend = agents.setup.pars['m_wage_trend'][t]
             
            #Arrays for preparation
            is_state = (state[:,t]==ist)
            ind = np.where(is_state)[0]
             
            if not np.any(is_state): continue
         
            zf,zm,psi=agents.setup.all_indices(t,iexo[ind,t])[1:4]
             
            #Relationship over time
            relt[ist,t]=np.sum(is_state)
             
            #Assets over time  
            ass_rel[ist,t]=np.mean(assets_t[ind,t])
             
            #Income over time
            if sname=="Female, single":
                inc_rel[ist,t]=np.mean(np.exp(agents.setup.exogrid.zf_t[t][zf]  + ftrend ))
                 
            elif sname=="Male, single":
                 inc_rel[ist,t]=np.mean(np.exp(agents.setup.exogrid.zf_t[t][zm] + mtrend))
                 
            elif sname=="Couple, C" or sname=="Couple, M":
                 inc_rel[ist,t]=np.mean(np.exp(agents.setup.exogrid.zf_t[t][zf] + ftrend)+np.exp(agents.setup.exogrid.zf_t[t][zm] + mtrend))
     
            else:
             
               print('Error: No relationship chosen')
                
    mdl.moments['share single'] = relt[0,:]/N
    mdl.moments['share mar'] = relt[2,:]/N
    mdl.moments['share coh'] = relt[3,:]/N
               
 
    nspells = (state[:,1:]!=state[:,:-1]).astype(np.int).sum(axis=1).max() + 1
     
    state_beg = -1*np.ones((N,nspells),dtype=np.int8)
    time_beg = -1*np.ones((N,nspells),dtype=np.bool)
    did_end = np.zeros((N,nspells),dtype=np.bool)
    state_end = -1*np.ones((N,nspells),dtype=np.int8)
    time_end = -1*np.ones((N,nspells),dtype=np.bool)
    sp_length = -1*np.ones((N,nspells),dtype=np.int16)
    is_spell = np.zeros((N,nspells),dtype=np.bool)
     
    state_beg[:,0] = 0 # THIS ASSUMES EVERYONE STARTS AS SINGLE
    time_beg[:,0] = 0
    sp_length[:,0] = 1
    is_spell[:,0] = True
    ispell = np.zeros((N,),dtype=np.int8)
     
    for t in range(1,agents.setup.pars['T']):
        ichange = (state[:,t-1] != state[:,t])
        sp_length[~ichange,ispell[~ichange]] += 1
         
        if not np.any(ichange): continue
         
        did_end[ichange,ispell[ichange]] = True
         
        is_spell[ichange,ispell[ichange]+1] = True
        sp_length[ichange,ispell[ichange]+1] = 1 # if change then 1 year right
        state_end[ichange,ispell[ichange]] = state[ichange,t]
        time_end[ichange,ispell[ichange]] = t-1
        state_beg[ichange,ispell[ichange]+1] = state[ichange,t] 
        time_beg[ichange,ispell[ichange]] = t
         
        ispell[ichange] = ispell[ichange]+1
         
         
    allspells_beg = state_beg[is_spell]
    allspells_len = sp_length[is_spell]
    allspells_end = state_end[is_spell] # may be -1 if not ended
     
    # If the spell did not end mark it as ended with the state at its start
    allspells_end[allspells_end==-1] = allspells_beg[allspells_end==-1]
     
    spells = np.stack((allspells_beg,allspells_len,allspells_end),axis=1)
     
     
     
     
    
    # TO FABIO: I did not edit things anything past this point.
    # Please look why my spells are different than yours.
     
     
    #Now divide spells by relationship nature
     
    for ist,sname in enumerate(state_codes):
        s = sname.replace(',', '')
        s = s.replace(' ', '')
         
         
        is_state= (spells[:,0]==ist)
         
        if not is_state.any(): continue
             
        globals()['spells_t'+s]=spells[is_state,:]
        is_state= (globals()['spells_t'+s][:,1]!=0)
        globals()['spells_'+s]=globals()['spells_t'+s][is_state,:]
     
     
    ##################################
    # Construct the Hazard functions
    #################################
         
    #Hazard of Divorce
    hazd=list()
    lgh=len(spells_CoupleM[:,0])
    for t in range(agents.setup.pars['T']):
         
        cond=spells_CoupleM[:,1]==t+1
        temp=spells_CoupleM[cond,2]
        cond1=temp!=2
        temp1=temp[cond1]
        if lgh>0:
            haz1=len(temp1)/lgh
            lgh=lgh-len(temp)
        else:
            haz1=0.0
        hazd=[haz1]+hazd
         
    hazd.reverse()
    hazd=np.array(hazd).T
     
    #Hazard of Separation
    hazs=list()
    lgh=len(spells_CoupleC[:,0])
    for t in range(agents.setup.pars['T']):
         
        cond=spells_CoupleC[:,1]==t+1
        temp=spells_CoupleC[cond,2]
        cond1=temp==0
        temp1=temp[cond1]
        if lgh>0:
            haz1=len(temp1)/lgh
            lgh=lgh-len(temp)
        else:
            haz1=0.0
        hazs=[haz1]+hazs
         
    hazs.reverse()
    hazs=np.array(hazs).T
     
    #Hazard of Marriage (Cohabitation spells)
    hazm=list()
    lgh=len(spells_CoupleC[:,0])
    for t in range(agents.setup.pars['T']):
         
        cond=spells_CoupleC[:,1]==t+1
        temp=spells_CoupleC[cond,2]
        cond1=temp==2
        temp1=temp[cond1]
        if lgh>0:
            haz1=len(temp1)/lgh
            lgh=lgh-len(temp)
        else:
            haz1=0.0
        hazm=[haz1]+hazm
         
    hazm.reverse()
    hazm=np.array(hazm).T
     
     
    mdl.moments['hazard sep'] = hazs
    mdl.moments['hazard div'] = hazd
    mdl.moments['hazard mar'] = hazm
     
     
 
     
    #Singles: Marriage vs. cohabitation transition
    #spells_s=np.append(spells_Femalesingle,spells_Malesingle,axis=0)
    spells_s = spells_Femalesingle
    cond=spells_s[:,2]>1
    spells_sc=spells_s[cond,2]
    condm=spells_sc==2
    sharem=len(spells_sc[condm])/max(len(spells_sc),0.0001)
     
    if draw:
     
        #Print something useful for debug and rest
        print('The share of singles choosing marriage is {0:.2f}'.format(sharem))
        cond=(state<2)
        if assets_t[cond].size:
            print('The max level of assets for singles is {:.2f}, the grid upper bound is {:.2f}'.format(np.amax(assets_t[cond]),max(agents.setup.agrid_s)))
        cond=(state>1)
        if assets_t[cond].size:
            print('The max level of assets for couples is {:.2f}, the grid upper bound is {:.2f}'.format(np.amax(assets_t[cond]),max(agents.setup.agrid_c)))
         
        #Setup a file for the graphs
        pdf = matplotlib.backends.backend_pdf.PdfPages("moments_graphs.pdf")
         
        #################
        #Get data moments
        #################
         
        #Get Data Moments
        with open('moments.pkl', 'rb') as file:
            packed_data=pickle.load(file)
         
            #Unpack Moments (see data_moments.py to check if changes)
            #(hazm,hazs,hazd,mar,coh,fls_ratio,W)
            hazm_d=packed_data[0]
            hazs_d=packed_data[1]
            hazd_d=packed_data[2]
            mar_d=packed_data[3]
            coh_d=packed_data[4]
            fls_d=np.ones(1)*packed_data[5]
 
         
         
        #############################################
        # Hazard of Divorce
        #############################################
        fig = plt.figure()
        f1=fig.add_subplot(2,1,1)
         
        plt.plot(np.array(range(len(hazd_d))), hazd[0:len(hazd_d)],'b',linewidth=1.5, label='Hazard of Divorce - S')
        plt.plot(np.array(range(len(hazd_d))), hazd_d,'r', linestyle='--',linewidth=1.5, label='Hazard of Divorce - D')
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3),
                  fancybox=True, shadow=True, ncol=3, fontsize='x-small')
        #plt.legend(loc='upper left', shadow=True, fontsize='x-small')
        plt.xlabel('Duration')
        plt.ylabel('Hazard')
         
        #############################################
        # Hazard of Separation
        #############################################
        fig = plt.figure()
        f1=fig.add_subplot(2,1,1)
         
        plt.plot(np.array(range(len(hazs_d))), hazs[0:len(hazs_d)],'b',linewidth=1.5, label='Hazard of Separation - S')
        plt.plot(np.array(range(len(hazs_d))), hazs_d,'r', linestyle='--',linewidth=1.5, label='Hazard of Separation - D')
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3),
                  fancybox=True, shadow=True, ncol=3, fontsize='x-small')
        #plt.legend(loc='upper left', shadow=True, fontsize='x-small')
        plt.xlabel('Duration')
        plt.ylabel('Hazard')
         
        #############################################
        # Hazard of Marriage
        #############################################
        fig = plt.figure()
        f1=fig.add_subplot(2,1,1)
         
        plt.plot(np.array(range(len(hazm_d))), hazm[0:len(hazm_d)],'b',linewidth=1.5, label='Hazard of Marriage - S')
        plt.plot(np.array(range(len(hazm_d))), hazm_d,'r', linestyle='--',linewidth=1.5, label='Hazard of Marriage - D')
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3),
                  fancybox=True, shadow=True, ncol=3, fontsize='x-small')
        #plt.legend(loc='upper left', shadow=True, fontsize='x-small')
        plt.xlabel('Duration')
        plt.ylabel('Hazard')
                    
        ##########################################
        # Assets Over the Live Cycle
        ##########################################
        fig = plt.figure()
        f2=fig.add_subplot(2,1,1)
         
        for ist,sname in enumerate(state_codes):
            plt.plot(np.array(range(agents.setup.pars['T'])), ass_rel[ist,],color=print(ist/len(state_codes)),markersize=6, label=sname)
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3),
                  fancybox=True, shadow=True, ncol=len(state_codes), fontsize='x-small')
        #plt.legend(loc='upper left', shadow=True, fontsize='x-small')
        plt.xlabel('Time')
        plt.ylabel('Assets')
         
        ##########################################
        # Income Over the Live Cycle
        ##########################################
        fig = plt.figure()
        f3=fig.add_subplot(2,1,1)
         
        for ist,sname in enumerate(state_codes):
           
            plt.plot(np.array(range(agents.setup.pars['T'])), inc_rel[ist,],color=print(ist/len(state_codes)),markersize=6, label=sname)
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3),
                  fancybox=True, shadow=True, ncol=len(state_codes), fontsize='x-small')
        plt.xlabel('Time')
        plt.ylabel('Income')
                 
                 
        ##########################################
        # Relationship Over the Live Cycle
        ##########################################      
        fig = plt.figure()
        f4=fig.add_subplot(2,1,1)
        for ist,sname in enumerate(state_codes):
            plt.plot([],[],color=print(ist/len(state_codes)), label=sname)
        plt.stackplot(np.array(range(agents.setup.pars['T'])),relt[0,]/N,relt[1,]/N,relt[2,]/N,relt[3,]/N,
                      colors = ['b','y','g','r'])           
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3),
                  fancybox=True, shadow=True, ncol=len(state_codes), fontsize='x-small')
        plt.xlabel('Time')
        plt.ylabel('Share')
         
        ##########################################
        # Relationship and Data
        ##########################################      
        fig = plt.figure()
        f4=fig.add_subplot(2,1,1)
        plt.plot(np.array(range(len(mar_d))), mar_d,'g',linewidth=1.5, label='Share Married - D')
        plt.plot(np.array(range(len(mar_d))), relt[2,0:len(mar_d)]/N,'g',linestyle='--',linewidth=1.5, label='Share Married - S')
        plt.plot(np.array(range(len(coh_d))), coh_d,'r',linewidth=1.5, label='Share Cohabiting - D')
        plt.plot(np.array(range(len(coh_d))), relt[3,0:len(coh_d)]/N,'r',linestyle='--',linewidth=1.5, label='Share Cohabiting - S')
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3),
                  fancybox=True, shadow=True, ncol=len(state_codes), fontsize='x-small')
        plt.xlabel('Time')
        plt.ylabel('Share')
         
        ##########################################
        # FLS Over the Live Cycle
        ##########################################      
        fig = plt.figure()
        f5=fig.add_subplot(2,1,1)
 
        plt.plot(np.array(range(agents.setup.pars['T'])), flsm,color='r', label='Marriage')
        plt.plot(np.array(range(agents.setup.pars['T'])), flsc,color='k', label='Cohabitation')         
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3),
                  fancybox=True, shadow=True, ncol=len(state_codes), fontsize='x-small')
        plt.xlabel('Time')
        plt.ylabel('FLS')
 
        ##########################################
        # Put graphs together
        ##########################################
        #show()
        for fig in range(1, plt.gcf().number + 1): ## will open an empty extra figure :(
            pdf.savefig( fig )
        
        pdf.close()
        matplotlib.pyplot.close("all")
         
        
