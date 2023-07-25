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
gen rel=0
replace rel=1 if (NUMCOHAB>=1 & NUMCOHAB<99) | (M95>=1 & M95<99)

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
gen month= date_c-birth if coh==1
replace month=date_m-birth if mar==1

gen clength=ENDCOH1-date_c if coh==1




*Sometimes cohabitation before marriage
*but cant remember the date_rele:correct this
*using date_rela at marriage as reference date_rele
*replace mar=0 if M113==1
*replace coh=1 if M113==1

*Sometimes they say that they cohabited with
*someone else other than husband before
*marriage. This creates some incongruences (sleep
*with another guy after marriage but before leaving together?)
*replace mar=0 if M116==1
*replace coh=1 if M116==1

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
replace unil=1935 if state=="Alaska" 
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
replace unil=1971 if state=="Idaho" 
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
replace unil=2009 if state=="Missouri" 
replace unil=1973 if state=="Montana" 
replace unil=1972 if state=="Nebraska" 
replace unil=1967 if state=="Nevada" 
replace unil=1971 if state=="New Hampshire" 
replace unil=2007 if state=="New Jersey" 
replace unil=1933 if state=="New Mexico" 
replace unil=2010 if state=="New York" 
replace unil=0 if state=="North Carolina" 
replace unil=1971 if state=="North Dakota" 
replace unil=1992 if state=="Ohio" 
replace unil=1953 if state=="Oklahoma" 
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


*Not interested in post 1988
replace unil=0 if unil>1988
replace eq=0 if unil>1988

*replace unil=unil+1 if unil!=0
gen unill=unil

*recode to make it comparable with date_reles
replace unil=(unil-1900)*12 if unil!=0
gen eqd=eq
replace eq=(eq-1900)*12 if eq!=0

*****************************************************







********************************************************************
*C) Final Adjustment
********************************************************************
*Round Some variables
replace birth = round(M485M/12+1900)
replace age_rel = int(age_rel)
gen date_relo = date_rel
replace date_rel=int(date_rel/12+1900)

replace M499A=100 if M499A>=56

gen native=1
replace native=0 if M499A>=56

replace state="Abroad" if M499A>=56






********************************************************************
*D) Sample Selection
********************************************************************

drop if native==0

*Dont consider cohabitaitons under one month
drop if clength<=0 & mar==0 & rel==1

*Keep if not too far away from policy change
*keep if (((date_rel-unill<=20) & (date_rel-unill>=-20)) | unill==0 )
*keep if (((date_rel-unill>=-20)) | unill==0 )


*What to do with age at first relationship<16?
*Many people fall here, check!
*keep if (birth>=1935 & birth<1988)
keep if date_rel>=1950 


gen agem=floor(M2BP01/10)
********************************************************************
*E) Controls
********************************************************************



********************************************************

*College
gen coll=0
replace coll=1 if COMPLED>=16

*Gender
gen female=1
replace female=0 if M2DP01==1

gen keep=0
replace keep=1 if M499B==1


*Numerical state
encode state,gen(state1)

*For children: dont intend to have+have
g nointc=0
replace nointc=1 if S88==2
g nch=M204



*Etnicity
gen etn=0
replace etn=1 if M484==1
replace etn=2 if M484==2
replace etn=3 if M484>=3



	
********************************************************************************
*5) Period equals one year
********************************************************************************
scalar myear=35
gen age_rel_old=age_rel
drop if age_rel<15
replace age_rel=myear if rel==0
replace age_rel=M2BP01  if  M2BP01<myear & rel==0
replace age_rel=myear if age_rel>myear & rel==1



gen dur=age_rel-15+1
gen id=_n
expand dur, generate(newvar)
bysort id: ge t = _n
sort id t


*Generate the end variable
sort id t
gen mar2=0
bysort id:replace mar2=1 if rel==1 & _n==_N  & mar==1 & age_rel_old<=myear

gen coh2=0
bysort id:replace coh2=1 if rel==1 & _n==_N  & mar==0 & age_rel_old<=myear 


*controls
*Month at which relationship started 
gen monthh=date_relo-(date_rel-1900)*12


*Time for unilateral divorce
gen unil2=unill
replace unil2=. if unill==0


*"Per" tells me how long one period is in years. If per was only one, the
*event study coefficient would be too jumpy
global per=1


gen cyear=birth+15+t

replace cyear=floor(cyear/$per)
replace date_rel=floor(date_rel/$per)
replace unil2=floor(unil2/$per)


********************************************************************************
********************************************************************************
***Actual analysis of divorce laws
********************************************************************************
********************************************************************************


********************************************************************************
*1) Basic regression
********************************************************************************

*Create controls
global controls      monthh state1 cyear coll etn t agem
global controls_weak state1 cyear

eststo clear
eststo:did_imputation mar2 id  cyear unil2 [weight=SAMWT] ,fe($controls_weak )      autosample cluster(state1)  maxit(200)
       su mar2 [weight=SAMWT]  if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	

eststo:did_imputation mar2 id  cyear unil2 [weight=SAMWT] ,fe($controls )      autosample cluster(state1)  maxit(200)
       su mar2 [weight=SAMWT]  if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	

	   
eststo:did_imputation coh2 id  cyear unil2 [weight=SAMWT] ,fe($controls_weak )      autosample cluster(state1) maxit(200)
       su coh2 [weight=SAMWT]  if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	

	   
eststo:did_imputation coh2 id  cyear unil2 [weight=SAMWT] ,fe($controls )      autosample cluster(state1) maxit(200)
       su coh2 [weight=SAMWT]  if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	





*Table with results
  #delimit ;
  esttab using "C:\Users\Fabio\Dropbox\py_mar_3\Analysis\output\reg_rel_smc.tex", label replace 
  b(4) se(4) star(* 0.10 ** 0.05 *** 0.01) stats(mD N, labels("Dependent variable mean" "Observations") fmt(%9.0f) ) keep(tau) varlabels(tau "Unilateral Divorce") eqlabels(None) nonotes depvars nonum nomtitles
  
  prehead(\begin{table}[H]\centering 
				
				\caption{The average effect of unilateral divorce on the choice of cohabitation vs. marriage vs. staying single for single individuals}  
				\label{tab:tabreld}
				\resizebox{\textwidth}{!}{
				\begin{tabular}{@{\extracolsep{4pt}}lcccc} \toprule
				 &\multicolumn{2}{c}{\textit{Dependent variable: marry (1), stay single (0) }}&\multicolumn{2}{c}{\textit{Dependent variable: cohabit (1), stay single (0)}}\\)	 
			 posthead(\cline{2-3}\cline{4-5})
             prefoot("State FE     & Yes & Yes& Yes & Yes\\
                      Year FE     & Yes & Yes& Yes & Yes\\
	                  Additional controls    & No & Yes& No& Yes\\")   
			postfoot(
				\bottomrule \noalign{\smallskip}
				\end{tabular}
				}
				\begin{minipage}{\textwidth}
				\scriptsize\smallskip  "Notes: This table reports the average effect of unilateral divorce on the choice of cohabitation vs. marriage vs. staying single using the methodology of \cite{borusyak2021} and data from the the \textit{relationships sample}. The unit of observation is one person-month where the respondent was single. Depending on the spefication, we study the risk of marriage (brekup) of singleness spells, where the event of breakup (marriage) is treated as censoring. The dummy variable \textit{Unilateral Divorce} takes value 1 if unilateral divorce was in place in state \textit{s} and year \textit{t} and 0 otherwise. The additional controls include dummies for ethnicity, age, education, and duration of the singleness spell. Standard errors are clustered at the state level. Coefficients that are significantly different from zero are denoted by *10\%, **5\%, and ***1\%."   \\ 
				\end{minipage}
				\end{table}
				)
				  ;
	#delimit cr


