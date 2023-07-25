*****************************************************
**NSFG88-construct data for analyzing cohabiting behavior
*****************************************************
clear all
use "C:\Users\Fabio\Dropbox\JMP\empirical analysis\NSFG\88\DS0001\09473-0001-Data.dta" 

*keep if ever cohabited
keep if COHEVER==1

*CONSTRUCT COHABITATION MEASURES

*Year cohabitation started
replace COHAB1=. if COHAB1==99797 | COHAB1==99898 | COHAB1==99999
replace COHAB1=COHAB1-90000 if COHAB1>90000 & COHAB1!=.
gen ysta=COHAB1
gen clength=COHABINT

*generate outcome
gen sep=0
replace sep=1 if COHOUT==4

gen mar=0
replace mar=1 if COHOUT==2 | COHOUT==3



********************************************************************
*generate a variable for indicating the date of unilateral divorce
*introduction and of equitable distribution(52 states)
********************************************************************
decode A_4,gen(state)

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

gen unill=unil

*recode to make it comparable with date_reles
replace unil=(unil-1900)*12 if unil!=0
gen eqd=eq
replace eq=(eq-1900)*12 if eq!=0


********************************************************************
*C) Additional Demographic Variables
********************************************************************

*College
gen coll=0
replace coll=1 if EDUCAT>=16


********************************************************************
*D) Final Adjustment
********************************************************************
*Round Some variables
gen date_relo=ysta
gen date_rel = int(ysta/12+1900)




gen native=1
replace native=0 if A_4>=56

replace state="Abroad" if A_4>=56

*RENAME VARIABELS FOR MAKING THEM COMPARABLE
gen keep=0
replace keep=1 if F_3==996

gen id=_n


****************************************************
*****************************************************
**Preparte to Expand the data***

gen clength1=clength+1
expand clength1, generate(newvar)
bysort id: ge t = _n
bysort id: ge faill = sep == 1 & _n==_N 
sort id t

*Variable for cohabitation end
ge fail2=0
sort id t
bysort id:replace fail2=1 if mar==1 & _n==_N
bysort id:replace fail2=2 if sep==1 & _n==_N
*****************************************************
*****************************************************

*VARIABLES FOR POLICY CHANGE****************
gen eqdd=0
replace eqdd=1 if date_rel>=eqd & eqd!=0

gen unid=0
replace unid=1 if date_rel>=unil & unil!=0

gen tit=0
replace tit=1 if date_rel<eqd 
replace tit=1 if eqd==0 & state=="Wisconsin" & date_rel<1986

gen com=0
replace com=1 if eqd==0
replace com=0 if eqd==0 & state=="Wisconsin" & date_rel<1986


*****************************************************

********************************************************************
*F) Sample Selection
********************************************************************


keep if native==1
/*
keep if rel<=3
keep if age_rel>=15
keep if age_rel<=65
keep if date_rel>=1950
keep if M2BP01<=70
keep if native==1
*/
drop if clength1<=0
keep if A_5<=75
keep if date_rel>=1950


*********************************************************************************
***G) Create controls
*********************************************************************************

*Numerical state
encode state,gen(state1)

*Risk of breakup
gen failc=0
replace failc=1 if fail2==2

*Risk of marriage
gen failm=0
replace failm=1 if fail2==1

*Risk of ending cohabitation
gen fail=0
replace fail=1 if failc+failm==1

*Time for unilateral divorce
gen unil2=unill
replace unil2=. if unill==0

*Dummy for unilateral divorce for couples that started before ud created
gen bau=0
replace bau=1 if date_rel<unil2 & date_rel+t/12>=unil2 & unil2!=.


gen agem=floor(A_5/10)

*Etnicity
gen etn=3
replace etn=1 if F9C==3
replace etn=2 if F9D==4


*Month at which relationship started 
gen monthh=date_relo-(date_rel-1900)*12

*Create dummy for current month
gen year=int((date_relo+t)/12+1900)
gen monthhh=date_relo+t-(year-1900)*12

*Duration of relatinoship per year
gen dur=floor(t/12)
replace dur=10 if dur>=10


*"Per" tells me how long one period is in years. If per was only one, the
*event study coefficient would be too jumpy
global per=1

gen date_rel_year=date_rel
gen unil2_year=unil2

replace date_rel=floor(date_rel/$per)
replace unil2=floor(unil2/$per)

*Dummy for unilateral divorce
gen unidd=0
replace unidd=1 if date_rel>=unil2 & unil2!=.

*interaction terms
gen int1=unidd*com
gen int2=unidd*tit
gen int3=unidd*eqdd


********************************************************************************
********************************************************************************
***Actual analysis of divorce laws
********************************************************************************
********************************************************************************

*******************************************************************************
*Note that some observation need to be dropped in did_estimation because
*some states are ALWAYS treated...
*******************************************************************************

*Create controls
global controls      bau dur monthh monthhh   state1 date_rel  coll etn agem
global controls_weak state1 date_rel 
global controls_sm      monthh state1 date_rel  coll etn agem
global controls_sm_weak state1 date_rel

eststo clear
eststo:did_imputation sep id  date_rel unil2 [weight=W5] if t==1 & sep+mar==1,fe($controls_sm_weak ) autosample cluster(state1)  
       su sep if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
eststo:did_imputation sep id  date_rel unil2 [weight=W5] if t==1 & sep+mar==1,fe($controls_sm )      autosample cluster(state1)  
       su sep if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
eststo:did_imputation sep id  date_rel unil2 [weight=W5] if t==1 & sep+mar==1 & keep==1,fe($controls_sm_weak ) autosample cluster(state1)  
       su sep if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
eststo:did_imputation sep id  date_rel unil2 [weight=W5] if t==1 & sep+mar==1 & keep==1,fe($controls_sm )      autosample cluster(state1)  
       su sep if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
	   
	   
	   
*Table with results
  #delimit ;
  esttab using "C:\Users\Fabio\Dropbox\py_mar_3\Analysis\output\reg_sep_nsfg.tex", label replace 
    b(4) se(4) star(* 0.10 ** 0.05 *** 0.01) stats(mD N, labels("Dependent variable mean" "N") fmt(%9.0f) ) keep(tau) varlabels(tau "Unilateral Divorce") eqlabels(None) nonotes depvars nomtitles nonum
  
  prehead(\begin{table}[H]\centering 
				\scriptsize
				\label{tab:tabsepnsfg}
				\caption{Likelihood that cohabitation ends in breakup vs. marriage. Unit of observation: period where cohabitation ends}  
				\resizebox{0.8\textwidth}{!}{
				\begin{tabular}{l*{@M}{c}} \toprule
				\textit{Dependent variable:}&\multicolumn{4}{c}{Marry (0), Breakup (1)}\\
				\textit{Sample:}& All &  All & NoMovers & NoMovers \\)	 	 
			 posthead(\hline)
             prefoot("\midrule state FE     & Yes & Yes& Yes & Yes\\
                      Period cohabitation began FE     & Yes & Yes& Yes & Yes\\
	                  Additional controls    & No & Yes& No& Yes\\")   
			postfoot(
				\bottomrule \noalign{\smallskip}
				\end{tabular}
				}
				\begin{minipage}{\textwidth}
				\scriptsize\smallskip Notes: "Standard errors are clustered at the state level. Coefficients that are significantly different from zero are denoted by *10\%, **5\% and ***1\%. Dataset: National Survey of Family Growth 1988 wave. The NoMovers sample includes only respondents who always resided in their current state of residence, while the All sample does not impose such a restriction. "  \\ 
				\end{minipage}
				\end{table}
				)
				  ;
	#delimit cr

	
*************************
**2 Cohabitation Duration
*************************
*Regressions
eststo clear
eststo:did_imputation fail id  date_rel unil2 [weight=W5],fe($controls_weak ) autosample cluster(state1)  
       su fail if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
	   
eststo:did_imputation fail id  date_rel unil2 [weight=W5],fe($controls )      autosample cluster(state1) 
       su fail if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	

eststo:did_imputation fail id  date_rel unil2 [weight=W5] if keep==1,fe($controls_weak ) autosample cluster(state1)
       su fail if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
	
eststo:did_imputation fail id  date_rel unil2 [weight=W5] if keep==1,fe($controls )      autosample cluster(state1) 
       su fail if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	



*Table with results
  #delimit ;
  esttab using "C:\Users\Fabio\Dropbox\py_mar_3\Analysis\output\reg_dur_nsfg.tex", label replace 
  b(4) se(4) star(* 0.10 ** 0.05 *** 0.01) stats(mD N, labels("Dependent variable mean" "N") fmt(%9.0f) ) keep(tau) varlabels(tau "Unilateral Divorce") eqlabels(None) nonotes depvars nomtitles nonum
  
  prehead(\begin{table}[H]\centering 
				\scriptsize
				\label{tab:tabdurnsfg}
				\caption{Likelihood that cohabitation spell ends. Unit of observation: cohabitation-month}  
				\resizebox{0.8\textwidth}{!}{
				\begin{tabular}{l*{@M}{c}} \toprule
				\textit{Dependent variable:}&\multicolumn{4}{c}{Keep cohabiting (0), cohabitation ends (1)}\\
				\textit{Sample:}& All &  All & NoMovers & NoMovers \\)	 
			 posthead(\hline)
             prefoot("\midrule state FE     & Yes & Yes& Yes & Yes\\
                      Period cohabitation began FE     & Yes & Yes& Yes & Yes\\
	                  Additional controls    & No & Yes& No& Yes\\")   
			postfoot(
				\bottomrule \noalign{\smallskip}
				\end{tabular}
				}
				\begin{minipage}{\textwidth}
				\scriptsize\smallskip Notes: "Standard errors are clustered at the state level. Coefficients that are significantly different from zero are denoted by *10\%, **5\% and ***1\%. Dataset: National Survey of Family Growth 1988 wave.  The NoMovers sample includes only respondents who always resided in their current state of residence, while the All sample does not impose such a restriction."  \\ 
				\end{minipage}
				\end{table}
				)
				  ;
	#delimit cr