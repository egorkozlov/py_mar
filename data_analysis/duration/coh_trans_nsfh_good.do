clear  all

import delimited using "C:\Users\Fabio\Dropbox\py_mar_3\histo.csv"



*drop if eq>1

gen res=numunion+numcohmr

gen inta=idatmm+(idatyy-1900)*12  


keep if (numunion-nummar>0) |  (numcohmr>0)  
drop if numcohmr==.
 
gen num=1


*keep and reshape
keep res  howbeg* begdat* enddat* howend* mardat* widdat* inta unil state eq birth

forvalues i=1(1)9{

rename howbeg0`i' howbeg`i'
rename enddat0`i' enddat`i'
rename howend0`i' howend`i'
rename mardat0`i' mardat`i'
rename begdat0`i' begdat`i'
rename widdat0`i' widdat`i'

}
* howend* cohmar* begdat* begflg* enddat* endflg* mardat* marflg* sepdat* sepflg* divdat* divflg* widdat* widflg* cohdat* cohflg*

drop howbeg9 howbeg10
drop enddat9 enddat10
drop howend9 howend10
drop mardat9 mardat10
drop begdat9 begdat10
drop widdat9 widdat10
gen id=_n

reshape long howbeg enddat howend widdat mardat begdat, i(id) j(rel)
*howend cohmar begdat begflg enddat endflg mardat marflg sepdat sepflg divdat divflg widdat widflg cohdat cohflg, i(id) j(rel)
keep if howbeg=="coh"

gen fine="censored"
replace fine="sep" if howend=="sep"
replace fine="mar" if howend=="div"
replace fine="mar" if howend=="intact" & mardat>1

gen end=-1
replace end=endd if fine=="sep"
replace end=mardat if fine=="mar"
replace end=widdat if howend=="wid"
replace end=inta if fine=="censored"


gen dur=end-begdat




gen dury=.

egen duryt = cut(dur), at(0,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,289,301,313,325,337,349,361,373,385,397)
replace dury=1 if duryt==0
replace dury=2 if duryt==13
replace dury=3 if duryt==25
replace dury=4 if duryt==37
replace dury=5 if duryt==49
replace dury=6 if duryt==61
replace dury=7 if duryt==73
replace dury=8 if duryt==85
replace dury=9 if duryt==97
replace dury=10 if duryt==109
replace dury=11 if duryt==121
replace dury=12 if duryt==133
replace dury=13 if duryt==145
replace dury=14 if duryt==157
replace dury=15 if duryt==169
replace dury=16 if duryt==181

replace dury=17 if duryt==193
replace dury=18 if duryt==205
replace dury=19 if duryt==217
replace dury=20 if duryt==229
replace dury=21 if duryt==241
replace dury=22 if duryt==253
replace dury=23 if duryt==265
replace dury=24 if duryt==277
replace dury=25 if duryt==289
replace dury=26 if duryt==301
replace dury=27 if duryt==313
replace dury=28 if duryt==325
replace dury=29 if duryt==337
replace dury=30 if duryt==349
replace dury=31 if duryt==361
replace dury=32 if duryt==373
replace dury=33 if duryt==385


drop dury
gen dury=round(dur/12,1)

*variable for unidivorce

gen begy = begdat/12+1900
replace begy=round(begy,1)
gen unid=0
replace unid=1 if begy>=unil

keep if dur>0 & dur<=2000
gen age=begy-birth
*keep if age>=14
*keep if birth>=1930
drop if begy>=2004

gen IDN1=id*1000+rel 


expand dur

encode state,gen(sta)

bysort id rel :gen t=_n
gen sepa=0
bysort id rel:replace sepa=1 if fin=="sep" & t==dur


gen mara=0
bysort id rel:replace mara=1 if fin=="mar" & t==dur

gen eqd=0
replace eqd=1 if begy>=eq & eq!=0


gen tit=0
replace tit=1 if begy<eq 


*Regression
reg sepa i.sta i.begy unid   if tit==0 & eqd==0,cluster(sta)


gen com=0
replace com=1 if eqd==0 & tit==0
gen int1=unid*com
gen int2=unid*eqd
gen int3=unid*tit

reg sepa i.sta i.begy  int1 int2 int3 eqd    if  begy<=1987,cluster(sta)


