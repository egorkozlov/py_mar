*****************************************************
**Build first Cohabitation Dataset
*****************************************************

*****************************************************
**1 NSFH
*****************************************************
clear  all
use "C:\Users\Fabio\Dropbox\JMP\empirical analysis\NSFH\Wave 1\ICPSR_06041\DS0001\06041-0001-Data.dta" 

*******************************************************
* A) Manupulation of Marital Hisory
******************************************************

*Keep if ever relationship
keep if (NUMCOHAB>=1 & NUMCOHAB<99) 

*Birth date
gen birth=M485M

*Year cohabitation started
gen date_c=BEGCOH1
replace date_c=. if BEGCOH1==0 |  BEGCOH1>2000


*Marriage dummy
gen mar=0
replace mar=1 if HOWEND1==2 | HOWEND1==3

*Separation Dummy
gen sep=0
replace sep=1 if HOWEND1==1


*date_rele at Interview
gen QDICMO=(MYEAR)*12+MMONTH

*Age in months at first relationship
drop if M485M>20000
gen monthc= date_c-birth 

*Cohabitation length
gen clength=ENDCOH1-date_c
replace clength=. if ENDCOH1==0 | ENDCOH1>2000


*date_rele relationship started
gen date_rel=BEGCOH1 

*Age at first Cohabitation
gen age_rel=(date_rel-birth)/12

*Id of the person
gen idp=_n



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


*recode to make it comparable with date_reles
replace unil=(unil-1900)*12+1 if unil!=0
replace eq=(eq-1900)*12+1 if eq!=0

*****************************************************



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



********************************************************************
*E) Adjust Time Fixed Effects
********************************************************************
bysort date_rel:egen marM=mean(mar)
gen date_rel_nsfh=date_rel

forval i=1/15{

replace date_rel_nsfh=date_rel_nsfh+1 if marM==1
drop marM
bysort date_rel_nsfh:egen marM=mean(mar)



}

saveold "C:\Users\Fabio\Dropbox\JMP\empirical analysis\NSFG\88\NSFHW1d.dta",replace



********************************************************************
*F) Expand the Data
********************************************************************
*replace length=length+1 if mar==1
*replace length=length+1 if sep==1
expand clength, generate(newvar)
bysort id: ge t = _n
bysort id: ge faill = sep == 1 & _n==_N 
sort id t

*Generate the end variable
ge fail2=0
sort id t
bysort id:replace fail2=1 if mar==1 & _n==_N
bysort id:replace fail2=2 if sep==1 & _n==_N
*****************************************************
*****************************************************

*VARIABLES FOR POLICY CHANGE****************
gen month=BEGCOH1+t-1
gen month1=BEGCOH1

gen eqd=0
replace eqd=1 if month>=eq & eq!=0

gen unid=0
replace unid=1 if month>=unil & unil!=0

gen unid1=0
replace unid1=1 if month1>=unil & unil!=0

gen tit=0
replace tit=1 if month<eq 

gen com=0
replace com=1 if eq==0
replace com=0 if eq==0 & state=="Wisconsin" & month<(1986-1900)*12+1

*interaction terms
gen int1=unid*com
gen int2=unid*tit
gen int3=unid*eqd

*******************************************************************


********************************************************************
*G) Sample Selection
********************************************************************

*drop if native==0

*Dont consider cohabitaitons under one month
drop if clength<=0

*Drop strange observations in NLSFH
drop if HOWEND1==3

*What to do with age at first relationship<16?
*Many people fall here, check!
*drop if age_rel<20 | age_rel>=.

*Adjust Weights
sum nsfh
gen Nnsfh=r(mean)
sum nsfh if nsfh==1
gen w_nsfh=r(sum)*r(N)
sum nsfg if nsfg==1
gen w_nsfg=r(sum)*r(N)
gen wgt_n=wgt*w_nsfh/(w_nsfh+w_nsfg)*(1/Nnsfh)
egen su=sum(wgt_n)
replace wgt_n=wgt_n/su


saveold "C:\Users\Fabio\Dropbox\JMP\empirical analysis\NSFG\88\NSFG88s.dta",replace

*********************************************************************************
*********************************************************************************
*********************************************************************************
