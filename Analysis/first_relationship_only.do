*****************************************************
**Build date_relaset of First Relationship
*****************************************************

*****************************************************
**1 NSFH
*****************************************************
clear all
use "C:\Users\Fabio\Dropbox\JMP\empirical analysis\NSFH\Wave 1\ICPSR_06041\DS0001\06041-0001-Data.dta" 

*******************************************************
* A) Manupulation of Marital Hisory
******************************************************

*Keep if ever relationship
keep if (NUMCOHAB>=1 & NUMCOHAB<99) | (M95>=1 & M95<99)

*Construct marriage and cohabitation measures
gen date_m=M96M if M96M<2000 

*Birth date
gen birth=M485M

*Year cohabitation started
gen date_c=BEGCOH1
replace date_c=. if BEGCOH1==0 |  BEGCOH1>2000

*Marriage dummy
gen mar=0
replace mar=1 if date_c>=date_m & date_c!=. & date_m!=. 
replace mar=1 if date_c==. 

*Cohabitation Dummy
gen coh=0
replace coh=1 if mar==0


*date_rele at Interview
gen QDICMO=(MYEAR)*12+MMONTH

*Age in months at first relationship
drop if M485M>20000
gen month= date_c-birth if coh==1
replace month=date_m-birth if mar==1

gen clength=ENDCOH1-date_c if coh==1




*Sometimes cohabitation before marriage
*but cant remember the date_rele:correct this
*using date_rela at marriage as reference date_rele
replace mar=0 if M113==1
replace coh=1 if M113==1

*Sometimes they say that they cohabited with
*someone else other than husband before
*marriage. This creates some incongruences (sleep
*with another guy after marriage but before leaving together?)
replace mar=0 if M116==1
replace coh=1 if M116==1

*date_rele relationship started
gen date_rel=BEGCOH1 if coh==1
replace date_rel=date_m if mar==1

*Age at first Marriage
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

*replace unil=unil+1 if unil!=0
gen unill=unil

*recode to make it comparable with date_reles
replace unil=(unil-1900)*12+1 if unil!=0
replace eq=(eq-1900)*12+1 if eq!=0

*****************************************************


*VARIABLES FOR POLICY CHANGE**************************

*Get the date_rele of singleness start
gen eqd=0
replace eqd=1 if date_rel>=eq & eq!=0

gen unid=0
replace unid=1 if date_rel>=unil & unil!=0

gen tit=0
replace tit=1 if date_rel<eq 

gen com=0
replace com=1 if eq==0
replace com=0 if eq==0 & state=="Wisconsin" & date_rel<(1986-1900)*12+1
replace com=0 if eq==0 & state=="Wisconsin" & month<(1986-1900)*12+1
********************************************************


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
gen birth_month=birth
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

*interaction terms
gen int1=unid*com
gen int2=unid*tit
gen int3=unid*eqd

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

saveold "C:\Users\Fabio\Dropbox\JMP\empirical analysis\NSFG\88\NSFHW1s.dta",replace


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
keep if COHEVER==1 | (F_13>=1 & F_13<=5)


*Construct marriage date
gen date_m=F19M1MO if F19M1MO!=99797 & F19M1MO!=99898 & F19M1MO!=99999 & F19M1MO!=0 
replace date_m=F15MO if BOX66==1 & F15MO!=99797 & F15MO!=99898 & F15MO!=99999 & F15MO!=0 & F15MO!=. & F19M1MO==0
replace date_m=F15MO if BOX66==1 & F15MO!=99797 & F15MO!=99898 & F15MO!=99999 & F15MO!=0 & F15MO!=. & F19M1MO==.
replace date_m=date_m-90000 if date_m>9000 

*Year cohabitation started
replace COHAB1=. if COHAB1==99797 | COHAB1==99898 | COHAB1==99999 
replace COHAB1=COHAB1-90000 if COHAB1>9000 
rename COHAB1 date_c


*Marriage variable
gen mar=0
replace mar=1 if date_c>=date_m & date_c!=. & date_m!=0 & date_c<2000 & date_m!=.  
replace mar=1 if date_c==. 
replace mar=1 if COHEVER==1 & COHSTAT==3 

*Cohabitation Variable
gen coh=0
replace coh=1 if  mar==0

*When relationship Started
gen date_rel=date_c if coh==1
replace date_rel=date_m if mar==1

*Birth
gen birth=A_3
*Age at first relationship
gen age_rel=(date_rel-birth)/12

*Date relationship starts
gen month= date_c-A_3 if coh==1
replace month=date_m-A_3 if mar==1
drop if month<0
gen length=month-12*20
drop if length<=0

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

*Generate the policy variables
gen eqd=0
replace eqd=1 if date_rel>=eq & eq!=0

gen unid=0
replace unid=1 if date_rel>=unil & unil!=0

gen tit=0
replace tit=1 if date_rel<eq 

gen com=0
replace com=1 if eq==0
replace com=0 if eq==0 & state=="Wisconsin" & date_rel<(1986-1900)*12+1

*interaction terms
gen int1=unid*com
gen int2=unid*tit
gen int3=unid*eqd

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

gen clength=COHABINT
********************************************************************
*D) Sample Selection
********************************************************************

*Drop if age at relationship crazy
drop if age_rel<=15

*RELACE IF COHABITATION SUPER SHORT**
drop if clength<1

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

********************************************************************
*E) Adjust Time Fixed Effects
********************************************************************
bysort date_rel:egen marM=mean(mar)
gen date_rel_nsfg=date_rel





forval i=1/15{

replace date_rel_nsfg=date_rel_nsfg+1 if marM==1
drop marM
bysort date_rel_nsfg:egen marM=mean(mar)



}


keep id   month birth_month unil eq age_rel  female  nsfg nsfh  mar coh  birth state unid date_rel  eqd tit com coll wgt int1 int2 int3 date_rel_nsfg clength native keep unill


***************************
*APPEND OTHER SURVEY NOW
**************************
append using "C:\Users\Fabio\Dropbox\JMP\empirical analysis\NSFG\88\NSFHW1s.dta"



********************************************************************
*F) Sample Selection
********************************************************************

*drop if native==0

*Dont consider cohabitaitons under one month
*drop if clength<=1 & mar==0 & reln==1

*What to do with age at first relationship<16?
*Many people fall here, check!
drop if age_rel<=15 | age_rel>=.
drop if clength<=0 & mar==0
drop if clength>300 & mar==0

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


*histograph begin
bysort stat:egen tot=mean(unid)
gen never=0
replace never=1 if tot==0
gen beg=date_rel


#delimit;
hist beg
         ,xlabel(1920(20)1988)
          xtitle(Year Relationship Began)
		  graphregion(color(white)) plotregion(color(white))
		   legend(row(1) label(1 "Simulated") ) ;
# delimit cr
graph export "C:\Users\Fabio\Dropbox\JMP\presentation\phd_apero_08_2019\hist_choice.pdf", as(pdf) replace


gen time=.
*replace time=1947.5 if beg>=1945 & beg<1950
replace time=1952.5 if beg>=1950 & beg<1955
replace time=1957.5 if beg>=1955 & beg<1960
replace time=1962.5 if beg>=1960 & beg<1965
replace time=1967.5 if beg>=1965 & beg<1970
replace time=1972.5 if beg>=1970 & beg<1975
replace time=1977.5 if beg>=1975 & beg<1980
replace time=1982.5 if beg>=1980 & beg<1985

bysort time never: egen mmar=mean(mar)

#delimit;
twoway (line mmar time if never==0 & mmar>0.5, color(gs8*.9) lwidth(0.45) yscale(r(0.6 1)) ) 
       (line mmar time if never==1 & mmar>0.5,color(maroon*.75) lwidth(0.45) yscale(r(0.6 1))),
	   ylabel(0.5(0.1)1)
graphregion(color(white)) plotregion(color(white)) xtitle(Relationship Began - Year) ytitle(Share Marriages) 
         legend(order(1 2) row(1) label(1 "Became Unilateral States") label(2 "Never Unilateral States") 
	label(3 "date_rela") );
# delimit cr
graph export "C:\Users\Fabio\Dropbox\JMP\presentation\phd_apero_08_2019\ev_choice.pdf", as(pdf) replace
*Reproduce result to get a graph


gen diffe=date_rel-unill






*binscatter res diffe1

*encode state,gen(sta)
*replace diffe1=diffe1+100
*reg mar i.date_rel i.sta ib90.diffe1 if diffe1!=.
*replace diffe1=diffe1-100



*Check of Parallel Trends
gen treated=1
replace treated=0 if unil==0
gen diffe_mod=diffe
replace diffe_mod=. if treated==0
reg  diffe i.date_rel i.birth if treated==1 & diffe<20 & diffe>-40
predict marh
replace diffe_mod=marh if treated==0

*binscatter mar diffe_mod if diffe1!=.,by(treated)

/*
gen diffe_mod1=.
replace diffe_mod1=0 if diffe_mod>=0 & diffe_mod<10
replace diffe_mod1=10 if diffe_mod>=10 & diffe_mod<20
*replace diffe_mod1=20 if diffe_mod>=20 & diffe_mod<30
replace diffe_mod1=-10 if diffe_mod>=-10 & diffe_mod<0
replace diffe_mod1=-20 if diffe_mod>=-20 & diffe_mod<-10
replace diffe_mod1=-30 if diffe_mod>=-30 & diffe_mod<-20
replace diffe_mod1=-40 if diffe_mod>=-40 & diffe_mod<-30
*replace diffe_mod1=-50 if diffe_mod>=-50 & diffe_mod<-40
*/

gen diffe_mod1=.
replace diffe_mod1=0 if diffe_mod>=0 & diffe_mod<5
replace diffe_mod1=5 if diffe_mod>=5 & diffe_mod<10
replace diffe_mod1=10 if diffe_mod>=10 & diffe_mod<15
*replace diffe_mod1=20 if diffe_mod>=20 & diffe_mod<30
replace diffe_mod1=-10 if diffe_mod>=-10 & diffe_mod<0
replace diffe_mod1=-20 if diffe_mod>=-20 & diffe_mod<-10
replace diffe_mod1=-30 if diffe_mod>=-30 & diffe_mod<-20
replace diffe_mod1=-40 if diffe_mod>=-40 & diffe_mod<-30
*replace diffe_mod1=-50 if diffe_mod>=-50 & diffe_mod<-40

sort diffe_mod1 treated

bysort diffe_mod1 treated:egen marm=mean(mar)

twoway (line marm diffe_mod1  if treated==0) (line marm diffe_mod1  if treated==1)

binscatter mar diffe_mod if diffe_mod<15 ,by(treated) xline(0) xtitle("Event Date") ytitle("Share Married") line(connect) legend(label(1 "Never Treated") label(2 "Treated")) savegraph("C:\Users\Fabio\Dropbox\JMP\empirical analysis\Analysis\ptrend.pdf") replace

reg mar i.date_rel if diffe_mod<3 & diffe_mod!=.
predict residual,resid

binscatter residual diffe_mod if diffe_mod<15 ,by(treated) 


