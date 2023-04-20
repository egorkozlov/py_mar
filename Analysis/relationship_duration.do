*****************************************************
*CONSTRUCT DURATION DATA FOR MULT LOGIT AND PROBIT
****************************************************
clear all
set maxvar 20000


*****************************************************
**1 NSFH
*****************************************************
clear all
use "C:\Users\Fabio\Dropbox\JMP\empirical analysis\NSFH\Wave 1\ICPSR_06041\DS0001\06041-0001-Data.dta" 


*******************************************************
* A) Manupulation of Marital Hisory
******************************************************




*Keep if ever relationship
keep if (NUMCOHAB==0) & (M95==0)

*Birth date
gen birth=M485M

*date_rele at Interview
gen QDICMO=(MYEAR)*12+MMONTH
gen date_rel=QDICMO

*Age in months at censoring
drop if M485M>20000
gen month= date_rel-birth

*Age at interview
gen age_rel=(date_rel-birth)/12

********************************************************************
*B) Generate Policy Variables
********************************************************************
decode M499A,gen(state)

*For policy variables
ge  mle=date_rel/12+1900
gen age=M2BP01

gen unil=0

replace unil=1971 if state=="Alabama" 
replace unil=1927 if state=="Alaska" 
replace unil=1973 if state=="Arizona" 
replace unil=0 if state=="Arkansas" 
replace unil=1970 if state=="California" 
replace unil=1972 if state=="Colorado" 
replace unil=1973 if state=="Connecticut" 
replace unil=1968 if state=="Delaware" 
replace unil=0 if state=="District of Columbia" 
replace unil=1971 if state=="Florida" 
replace unil=1973 if state=="Georgia" 
replace unil=1972 if state=="Hawaii" 
replace unil=1961 if state=="Idaho" 
replace unil=0 if state=="Illinois" 
replace unil=1973 if state=="Indiana" 
replace unil=1970 if state=="Iowa" 
replace unil=1969 if state=="Kansas" 
replace unil=1972 if state=="Kentucky" 
replace unil=0 if state=="Louisiana" 
replace unil=1973 if state=="Maine" 
replace unil=0 if state=="Maryland" 
replace unil=1975 if state=="Massachusetts" 
replace unil=1972 if state=="Michigan" 
replace unil=1974 if state=="Minnesota" 
replace unil=0 if state=="Mississippi" 
replace unil=0 if state=="Missouri" 
replace unil=1973 if state=="Montana" 
replace unil=1972 if state=="Nebraska" 
replace unil=1967 if state=="Nevada" 
replace unil=1971 if state=="New Hampshire" 
replace unil=0 if state=="New Jersey" 
replace unil=1927 if state=="New Mexico" 
replace unil=0 if state=="New York" 
replace unil=0 if state=="North Carolina" 
replace unil=1971 if state=="North Dakota" 
replace unil=1992 if state=="Ohio" 
replace unil=1927 if state=="Oklahoma" 
replace unil=1971 if state=="Oregon" 
replace unil=0 if state=="Pennsylvania" 
replace unil=1975 if state=="Rhode Island" 
replace unil=0 if state=="South Carolina" 
replace unil=1985 if state=="South Dakota" 
replace unil=0 if state=="Tennessee" 
replace unil=1970 if state=="Texas" 
replace unil=1987 if state=="Utah" 
replace unil=0 if state=="Vermont" 
replace unil=0 if state=="Virginia" 
replace state="Washington state" if state=="Washington"
replace unil=1973 if state=="Washington state" 
replace unil=1984 if state=="West Virginia" 
replace unil=1978 if state=="Wisconsin" 
replace unil=1977 if state=="Wyoming" 


gen eq=0

replace eq=1984 if state=="Alabama" 
replace eq=1927 if state=="Alaska" 
replace eq=0 if state=="Arizona" 
replace eq=1977 if state=="Arkansas" 
replace eq=0 if state=="California" 
replace eq=1972 if state=="Colorado" 
replace eq=1973 if state=="Connecticut" 
replace eq=1927 if state=="Delaware" 
replace eq=1977 if state=="District of Columbia" 
replace eq=1980 if state=="Florida" 
replace eq=1984 if state=="Georgia" 
replace eq=1927 if state=="Hawaii" 
replace eq=0 if state=="Idaho" 
replace eq=1977 if state=="Illinois" 
replace eq=1927 if state=="Indiana" 
replace eq=1927 if state=="Iowa" 
replace eq=1927 if state=="Kansas" 
replace eq=1976 if state=="Kentucky" 
replace eq=0 if state=="Louisiana" 
replace eq=1972 if state=="Maine" 
replace eq=1978 if state=="Maryland" 
replace eq=1974 if state=="Massachusetts" 
replace eq=1927 if state=="Michigan" 
replace eq=1927 if state=="Minnesota" 
replace eq=1989 if state=="Mississippi" 
replace eq=1977 if state=="Missouri" 
replace eq=1976 if state=="Montana" 
replace eq=1972 if state=="Nebraska" 
replace eq=0 if state=="Nevada" 
replace eq=1977 if state=="New Hampshire" 
replace eq=1974 if state=="New Jersey" 
replace eq=1967 if state=="New Mexico" 
replace eq=1980 if state=="New York" 
replace eq=1981 if state=="North Carolina" 
replace eq=1927 if state=="North Dakota" 
replace eq=1981 if state=="Ohio" 
replace eq=1975 if state=="Oklahoma" 
replace eq=1971 if state=="Oregon" 
replace eq=1980 if state=="Pennsylvania" 
replace eq=1981 if state=="Rhode Island" 
replace eq=1985 if state=="South Carolina" 
replace eq=1927 if state=="South Dakota" 
replace eq=1927 if state=="Tennessee" 
replace eq=0 if state=="Texas" 
replace eq=1927 if state=="Utah" 
replace eq=1927 if state=="Vermont" 
replace eq=1982 if state=="Virginia" 
replace eq=0 if state=="Washington state" 
replace eq=1985 if state=="West Virginia" 
replace eq=0 if state=="Wisconsin" 
replace eq=1967 if state=="Wyoming" 

*replace unil=unil+1 if unil!=0
gen unill=unil

*recode to make it comparable with date_reles
replace unil=(unil-1900)*12+1 if unil!=0
replace eq=(eq-1900)*12+1 if eq!=0


********************************************************************
*C) Additional Demographic Variables
********************************************************************

*College
gen coll=0
replace coll=1 if COMPLED>=16

*Gender
gen female=1
replace female=0 if M2DP01==1

********************************************************************
*D) Final Adjustment
********************************************************************
*Round Some variables
gen birth_month=M485M
replace birth = round(M485M/12+1900)
replace age_rel = round(age_rel)
replace date_rel = round(date_rel/12+1900)

replace M499A=100 if M499A>=56

gen native=1
replace native=0 if M499A>=56

replace state="Abroad" if M499A>=56

*RENAME VARIABELS FOR MAKING THEM COMPARABLE
gen nsfh=1
gen nsfg=0
gen keep=0
replace keep=1 if M499B==1
gen id=_n
replace id=id+1000000
rename SAMWT wgt



saveold "C:\Users\Fabio\Dropbox\JMP\empirical analysis\NSFG\88\NSFHW1sin.dta",replace

*****************************************************
**1 NSFG88
*****************************************************
clear all
set matsize 10000
use "C:\Users\Fabio\Dropbox\JMP\empirical analysis\NSFG\88\DS0001\09473-0001-Data.dta" 

*******************************************************
* A) Manupulation of Marital Hisory
******************************************************

*keep if in a relationship
keep if COHEVER==2 & (BOX64>=4)

*Birth
gen birth=A_3

*Date at interview
gen date_rel=QDICMO

*Age at first relationship
gen age_rel=(date_rel-birth)/12

*Date relationship starts
gen month= date_rel-A_3 

********************************************************************
*B) Generate Policy Variables
********************************************************************
decode A_4,gen(state)

gen unil=0

replace unil=1971 if state=="Alabama" 
replace unil=1927 if state=="Alaska" 
replace unil=1973 if state=="Arizona" 
replace unil=0 if state=="Arkansas" 
replace unil=1970 if state=="California" 
replace unil=1972 if state=="Colorado" 
replace unil=1973 if state=="Connecticut" 
replace unil=1968 if state=="Delaware" 
replace unil=0 if state=="District of Columbia" 
replace unil=1971 if state=="Florida" 
replace unil=1973 if state=="Georgia" 
replace unil=1972 if state=="Hawaii" 
replace unil=1961 if state=="Idaho" 
replace unil=0 if state=="Illinois" 
replace unil=1973 if state=="Indiana" 
replace unil=1970 if state=="Iowa" 
replace unil=1969 if state=="Kansas" 
replace unil=1972 if state=="Kentucky" 
replace unil=0 if state=="Louisiana" 
replace unil=1973 if state=="Maine" 
replace unil=0 if state=="Maryland" 
replace unil=1975 if state=="Massachusetts" 
replace unil=1972 if state=="Michigan" 
replace unil=1974 if state=="Minnesota" 
replace unil=0 if state=="Mississippi" 
replace unil=0 if state=="Missouri" 
replace unil=1973 if state=="Montana" 
replace unil=1972 if state=="Nebraska" 
replace unil=1967 if state=="Nevada" 
replace unil=1971 if state=="New Hampshire" 
replace unil=0 if state=="New Jersey" 
replace unil=1927 if state=="New Mexico" 
replace unil=0 if state=="New York" 
replace unil=0 if state=="North Carolina" 
replace unil=1971 if state=="North Dakota" 
replace unil=1992 if state=="Ohio" 
replace unil=1927 if state=="Oklahoma" 
replace unil=1971 if state=="Oregon" 
replace unil=0 if state=="Pennsylvania" 
replace unil=1975 if state=="Rhode Island" 
replace unil=0 if state=="South Carolina" 
replace unil=1985 if state=="South Dakota" 
replace unil=0 if state=="Tennessee" 
replace unil=1970 if state=="Texas" 
replace unil=1987 if state=="Utah" 
replace unil=0 if state=="Vermont" 
replace unil=0 if state=="Virginia" 
replace state="Washington state" if state=="Washington"
replace unil=1973 if state=="Washington state" 
replace unil=1984 if state=="West Virginia" 
replace unil=1978 if state=="Wisconsin" 
replace unil=1977 if state=="Wyoming" 


gen eq=0

replace eq=1984 if state=="Alabama" 
replace eq=1927 if state=="Alaska" 
replace eq=0 if state=="Arizona" 
replace eq=1977 if state=="Arkansas" 
replace eq=0 if state=="California" 
replace eq=1972 if state=="Colorado" 
replace eq=1973 if state=="Connecticut" 
replace eq=1927 if state=="Delaware" 
replace eq=1977 if state=="District of Columbia" 
replace eq=1980 if state=="Florida" 
replace eq=1984 if state=="Georgia" 
replace eq=1927 if state=="Hawaii" 
replace eq=0 if state=="Idaho" 
replace eq=1977 if state=="Illinois" 
replace eq=1927 if state=="Indiana" 
replace eq=1927 if state=="Iowa" 
replace eq=1927 if state=="Kansas" 
replace eq=1976 if state=="Kentucky" 
replace eq=0 if state=="Louisiana" 
replace eq=1972 if state=="Maine" 
replace eq=1978 if state=="Maryland" 
replace eq=1974 if state=="Massachusetts" 
replace eq=1927 if state=="Michigan" 
replace eq=1927 if state=="Minnesota" 
replace eq=1989 if state=="Mississippi" 
replace eq=1977 if state=="Missouri" 
replace eq=1976 if state=="Montana" 
replace eq=1972 if state=="Nebraska" 
replace eq=0 if state=="Nevada" 
replace eq=1977 if state=="New Hampshire" 
replace eq=1974 if state=="New Jersey" 
replace eq=1967 if state=="New Mexico" 
replace eq=1980 if state=="New York" 
replace eq=1981 if state=="North Carolina" 
replace eq=1927 if state=="North Dakota" 
replace eq=1981 if state=="Ohio" 
replace eq=1975 if state=="Oklahoma" 
replace eq=1971 if state=="Oregon" 
replace eq=1980 if state=="Pennsylvania" 
replace eq=1981 if state=="Rhode Island" 
replace eq=1985 if state=="South Carolina" 
replace eq=1927 if state=="South Dakota" 
replace eq=1927 if state=="Tennessee" 
replace eq=0 if state=="Texas" 
replace eq=1927 if state=="Utah" 
replace eq=1927 if state=="Vermont" 
replace eq=1982 if state=="Virginia" 
replace eq=0 if state=="Washington state" 
replace eq=1985 if state=="West Virginia" 
replace eq=0 if state=="Wisconsin" 
replace eq=1967 if state=="Wyoming" 
gen unill=unil

*replace unil=unil+1 if unil!=0

*recode to make it comparable with date_reles
replace unil=(unil-1900)*12+1 if unil!=0
replace eq=(eq-1900)*12+1 if eq!=0


********************************************************************
*C) Additional Demographic Variables
********************************************************************
*education
gen coll=0
replace coll=1 if EDUCAT>=16

gen native=1
replace native=0 if A_4>56

replace state="Abroad" if A_4>56

gen birth_month=birth

gen female=1


********************************************************************
*D) Sample Selection
********************************************************************


*NON AMERICA BORN
*keep if native==1

********************************************************************
*E) Finale Adjustments
********************************************************************
replace birth = round(A_3/12+1900)
replace age_rel = round(age_rel)
replace date_rel = round(date_rel/12+1900)
replace date_rel=round(birth+age_rel)
ge id=_n
gen nsfh=0
gen nsfg=1
gen reln=1

gen keep=0
replace keep=1 if F_3==996

gen wgt=W1

***************************
*APPEND OTHER SURVEY NOW
**************************
append using "C:\Users\Fabio\Dropbox\JMP\empirical analysis\NSFG\88\NSFHW1sin.dta"

g mar=0
g coh=0
global var id   state month birth_month unil eq age_rel  female  nsfg nsfh  mar coh  birth state  date_rel  coll wgt  native keep unill
keep $var

*Now put together relationship with relationship data
append using "C:\Users\Fabio\Dropbox\JMP\empirical analysis\NSFG\88\NSFG88s.dta",keep($var)


*******************************
*MAKE DATA OF THE DURATION Type
*******************************
gen length=month-15*12
drop if length<=0

gen idd=_n
*expand length
expand length
bysort idd: gen _t=_n
gen time=_t+15*12+birth_month
sort idd time

*******************************
*Generate the policy variables
*******************************
gen eqd=0
replace eqd=1 if time>=eq & eq!=0

gen unid=0
replace unid=1 if time>=unil & unil!=0

gen tit=0
replace tit=1 if time<eq 

gen com=0
replace com=1 if eq==0
replace com=0 if eq==0 & state=="Wisconsin" & time<(1986-1900)*12+1


*interaction terms
gen int1=unid*com
gen int2=unid*tit
gen int3=unid*eqd
gen int4=unid*(eqd+com)

**********************************************************************
*Keep only one year backwards to avoid having too many observations
*********************************************************************
gen itime=-(_t-length)

*Keep only years
g month_mod=mod(time,12)
replace month_mod=0 if month_mod<=0
gen timem=time-month_mod

*Make year
gen timey=int((time-month_mod)/12+1900)

g relt=0
replace relt=1 if coh==1 & itime==0
replace relt=2 if mar==1 & itime==0

bysort idd timey:egen rel=max(relt)

keep if month_mod==0
*bysort idd timey:keep if _n==1



drop id
rename idd id

gen id2=_n
encode state,gen(st)


sort id timey
bysort id :gen duration=_n
gen duration2=duration^2
gen duration3=duration^3
keep if timey>=1955

mlogit rel unid i.st  i.timey duration duration2 coll female                    ,baseoutcome(2) rrr
mlogit rel unid i.st  i.timey duration duration2 coll female if keep==1         ,baseoutcome(2) rrr
mlogit rel int4 int2 tit i.st  i.timey duration duration2 coll female           ,baseoutcome(2) rrr
mlogit rel int4 int2 tit i.st  i.timey duration duration2 coll female if keep==1,baseoutcome(2) rrr

/*
keep if birth>1947
quietly:mlogit rel unid i.st i.timey coll female duration,baseoutcome(2)
predict p
quietly:mlogit rel unid i.st i.timey coll female duration if keep==1,baseoutcome(2)
predict p1
quietly:mlogit rel unid i.st i.timey coll female duration  if nsfh==1,baseoutcome(2)
predict p2
quietly:mlogit rel unid i.st i.timey coll female duration if nsfh==0,baseoutcome(2)
predict p3

saveold "C:\Users\Fabio\Dropbox\JMP\empirical analysis\NSFG\88\NSFG88s2.dta",replace

/*
*drop if birth<1946
*drop if birth>1972 
*drop if duration>30
mlogit rel unid  i.st i.timey c.timey#st, baseoutcome(2)
predict p
mlogit rel unid  i.st i.timey if nsfh==0, baseoutcome(2)
predict p1
*drop if p>0.999
*drop if p1>0.999

/*


mlogit rel unid i.st i.timey, baseoutcome(2) cluster(st)
mlogit rel int2 int4 tit i.st i.timey, baseoutcome(2) cluster(st)

**********************************************************************
*Keep only one year backwards to avoid having too many observations
*********************************************************************
gen itime=-(_t-length)

*Keep only years
keep if mod(itime,12)==0

g rel=0
replace rel=1 if coh==1 & itime==0
replace rel=2 if mar==1 & itime==0

gen timey=round(time/12)

drop id
rename idd id

gen id2=_n
encode state,gen(st)

saveold "C:\Users\Fabio\Dropbox\JMP\empirical analysis\NSFG\88\NSFG88s2.dta",replace

