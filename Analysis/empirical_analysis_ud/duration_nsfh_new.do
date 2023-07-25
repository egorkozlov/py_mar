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




*keep and reshape
keep BEGCOH* ENDCOH* CONTIN* HOWEND*  birth  MYEAR MMONTH M485M M499A M2BP01 COMPLED M2DP01 M485M M499A M499B SAMWT NUMCOHAB M484 EDUCAT T58A  M488 M424B M424A


gen id2=_n

reshape long BEGCOH ENDCOH CONTIN HOWEND, i(id2) j(rel)

keep if rel<=NUMCOHAB



*Marriage dummy
gen mar=0
replace mar=1 if HOWEND==2 | HOWEND==3

*Separation Dummy
gen sep=0
replace sep=1 if HOWEND==1

*Cohabitation length
gen clength=ENDCOH-BEGCOH



*Age at first Cohabitation
gen age_rel=(BEGCOH-birth)/12




********************************************************************
*B) Generate Policy Variables
********************************************************************
decode M499A,gen(state)

*For policy variables
ge  mle=BEGCOH/12+1900
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
gen date_relo=BEGCOH
gen date_rel = int(BEGCOH/12+1900)

replace M499A=100 if M499A>=56

gen native=1
replace native=0 if M499A>=56

replace state="Abroad" if M499A>=56



*RENAME VARIABELS FOR MAKING THEM COMPARABLE
gen keep=0
replace keep=1 if M499B==1

gen id=_n



********************************************************************
*D) Expand the Data
********************************************************************
replace clength=clength+1
expand clength, generate(newvar)
bysort id: ge t = _n
bysort id: ge faill = sep == 1 & _n==_N 
gen idt=_n
sort id t


*Generate the end variable
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



*******************************************************************


********************************************************************
*E) Sample Selection
********************************************************************

*drop if native==0

*Dont consider cohabitaitons under one month
drop if  ENDCOH>=2000
drop if HOWEND>4
drop if BEGCOH==0 |  BEGCOH>2000


drop if clength<=0
keep if native==1
keep if date_rel>=1950
/*
keep if rel<=3
keep if age_rel>=15
keep if age_rel<=65
keep if date_rel>=1950
keep if M2BP01<=70
keep if native==1
*/

*keep if M2BP01<=75
*Drop strange observations in NLSFH
*drop if HOWEND==3

*Keep if not too far away from policy change
*keep if (((date_rel-unill<=20) & (date_rel-unill>=-20)) | unill==0 )
*keep if (((date_rel-unill>=-20)) | unill==0 )

*********************************************************************************
***F) Create controls
*********************************************************************************

*College
gen coll=0
replace coll=1 if COMPLED>=16

*Gender
gen female=1
replace female=0 if M2DP01==1


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

gen agem=floor(M2BP01/10)

*Etnicity
gen etn=0
replace etn=1 if M484==1
replace etn=2 if M484==2
replace etn=3 if M484>=3

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
global per=2

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
global controls      bau dur    state1 date_rel coll etn agem
global controls_weak state1 date_rel 
global controls_year      bau dur   state1 date_rel_year coll etn  agem
global controls_weak_year state1 date_rel_year
global controls_sm      state1 date_rel coll etn agem
global controls_sm_weak state1 date_rel
global controls_sm_year      state1 date_rel coll etn agem
global controls_sm_weak_year state1 date_rel
global nos state1!=51 | state1!=40 | state1!=35 | state1!=27 | state1!=28


**********************
**0 Summary statistics
*********************
gen censored=0
replace censored=1 if mar==0 & sep==0

label variable censored "End of Cohabitation spell: Breakup (0/1)"
label variable mar "End of Cohabitation spell: Marry (0/1)"
label variable sep "End of Cohabitation spell: Censored (0/1)"
label variable M2BP01 "Age of respondent"
label variable coll "College (0/1)"
label variable date_rel_year "Year the Cohabitation began"
label variable unidd "Unilateral divorce when Cohabitation began (0/1)"
label variable clength "Length of Cohabitation (months)"

est clear
estpost tabstat sep mar censored coll clength date_rel_year unidd M2BP01 [weight=SAMWT] if t==1, c(stat) stat(count mean median sd)
#delimit ;
esttab using "C:\Users\Fabio\Dropbox\py_mar_3\Analysis\output\sum_coh.tex", replace
 cells("count mean(fmt(%6.2fc)) p50 sd(fmt(%6.2fc))")   nonumber 
  nomtitle nonote noobs label booktabs 
  collabels("N" "Mean" "Median" "SD") ;
#delimit cr
	


**********************
**1 Cohabitation End
*********************

*Regressions
eststo clear
eststo:did_imputation sep id  date_rel_year unil2_year [weight=SAMWT] if t==1 & sep+mar==1,fe($controls_sm_weak_year ) autosample cluster(state1)  
       su sep if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
eststo:did_imputation sep id  date_rel_year unil2_year [weight=SAMWT] if t==1 & sep+mar==1,fe($controls_sm_year )      autosample cluster(state1)  
       su sep if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
eststo:did_imputation sep id  date_rel_year unil2_year [weight=SAMWT] if t==1 & sep+mar==1 & $nos,fe($controls_sm_weak_year ) autosample cluster(state1)   timecontrols(state1) 
       su sep if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
eststo:did_imputation sep id  date_rel_year unil2_year [weight=SAMWT] if t==1 & sep+mar==1  & $nos,fe($controls_sm_year )      autosample cluster(state1)   timecontrols(state1) 
       su sep if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	


*Table with results
  #delimit ;
  esttab using "C:\Users\Fabio\Dropbox\py_mar_3\Analysis\output\reg_sep.tex", label replace 
    b(4) se(4) star(* 0.10 ** 0.05 *** 0.01) stats(mD N, labels("Dependent variable mean" "Observations") fmt(%9.0f) ) keep(tau) varlabels(tau Unilateral Divorce}) eqlabels(None) nonotes depvars nomtitles nonum
  
  prehead(\begin{table}[H]\centering 
				\scriptsize
				\caption{The average effect of unilateral divorce on the probability that cohabitation ends in breakup vs. marriage. Unit of observation: end of a cohabitation spell}  
				\label{tab:tabsep}
				\resizebox{0.8\textwidth}{!}{
				\begin{tabular}{l*{@M}{c}} \toprule
			    &\multicolumn{4}{c}{\textit{Dependent variable: marry (0), breakup (1)}}\\)	 
			 posthead(\midrule)
             prefoot("State FE     & Yes & Yes& Yes & Yes\\
                      Year cohabitation began FE     & Yes & Yes& Yes & Yes\\
	                  Additional controls    & No & Yes& No& Yes\\
					  State-specific linear time trend    & No & No & Yes & Yes\\")   
			postfoot(
				\bottomrule \noalign{\smallskip}
				\end{tabular}
				}
				\begin{minipage}{\textwidth}
				\scriptsize\smallskip "Notes: This table reports the average effect of unilateral divorce on the probability of cohabitation ending in a breakup vs. marriage. The analysis follows the methodology outlined in \cite{borusyak2021} uses and the \textit{cohabitations sample}. The unit of observation is the end of a cohabitation spell (either marriage or breakup) which \textit{started} in state \textit{s} and year \textit{t}. The dummy variable \textit{Unilateral Divorce} takes value 1 if unilateral divorce was in effect in state \textit{s} and year \textit{t}. The additional controls include dummies for ethnicity, age, and education. The sample associated to regressions including State-specific linear time trend excludes five states due to multicollinearity. Standard errors are clustered at the state level. Coefficients that are significantly different from zero are denoted by *10\%, **5\%, and ***1\%."  \\ 
				\end{minipage}
				\end{table}
				)
				  ;
	#delimit cr
	
	


/*
*Event studies
did_imputation sep id  date_rel_year unil2_year [weight=SAMWT] if t==1 & sep+mar==1,            fe($controls_sm_weak ) autosample cluster(state1) maxit(200) allhorizons pre(8)
estimates store sep1
did_imputation sep id  date_rel_year unil2_year [weight=SAMWT] if t==1 & sep+mar==1,            fe($controls_sm )      autosample cluster(state1) maxit(200) allhorizons pre(8)
estimates store sep2
did_imputation sep id  date_rel_year unil2_year [weight=SAMWT] if t==1 & sep+mar==1   & $nos,           fe($controls_sm_weak ) autosample cluster(state1) maxit(200) allhorizons pre(8)
estimates store sep3
did_imputation sep id  date_rel_year unil2_year [weight=SAMWT] if t==1 & sep+mar==1   & $nos,           fe($controls_sm )      autosample cluster(state1) maxit(200) allhorizons pre(8)
estimates store sep4


#delimit ;
event_plot sep1 sep2 sep3 sep4, 
     plottype(scatter) ciplottype(rcap) 
	together perturb(-0.325(0.13)0.325) trimlead(7) trimlag(8) noautolegend 
	graph_opt(
		xtitle("Years since the event") ytitle("Average causal effect")  yscale(titlegap(*-35))
		xlabel(-7 "(-14,-13)" -6 "(-12,-11)" -5 "(-10,-9)" -4 "(-8,-7)" -3 "(-6,-5)" -2 "(-4,-3)" -1 "(-2,-1)"
		        0 "(0,1)"      1 "(2,3)"     2 "(4,5)"    3 "(6,7)"    4 "(8,9)"    5 "(10,11) "   6 "(10,11)"  7 "(12,13)" 8 "(14,15)", format(%18.0fc) labsize(6pt)) ylabel(-0.6(0.1)0.6) 
		legend(order(1 "Basic controls" 3 "Full controls" 5 "Basic controls  + State-specific linear time trends" 7 "Full controls + State-specific linear time trends") rows(4) region(style(none))) 
		xline(-0.5, lcolor(gs8) lpattern(dash)) yline(0, lcolor(gs8)) graphregion(color(white)) bgcolor(white) ylabel(, angle(horizontal)) 
	) 
	lag_opt1(msymbol(Oh)  color(cranberry))    lag_ci_opt1(color(cranberry)) 
	lag_opt2(msymbol(Sh) color(forest_green)) lag_ci_opt2(color(forest_green)) 	
	lag_opt3(msymbol(Dh) color(navy))         lag_ci_opt3(color(navy))
	lag_opt4(msymbol(Th) color(dkorange))     lag_ci_opt4(color(dkorange));
#delimit cr
graph export "C:\Users\Fabio\Dropbox\py_mar_3\Analysis\output\event_graph_sep.png",replace
*/
*************************
**2 Cohabitation Duration
*************************
*Regressions
eststo clear
eststo:did_imputation fail id  date_rel_year unil2_year [weight=SAMWT],fe($controls_weak_year ) autosample cluster(state1)  
       su fail if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
	   
eststo:did_imputation fail id  date_rel_year unil2_year [weight=SAMWT],fe($controls_year )      autosample cluster(state1) 
       su fail if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	

eststo:did_imputation fail id  date_rel_year unil2_year [weight=SAMWT],fe($controls_weak_year ) autosample cluster(state1) timecontrols(state1) 
       su fail if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
	
eststo:did_imputation fail id  date_rel_year unil2_year [weight=SAMWT],fe($controls_year )      autosample cluster(state1) timecontrols(state1)   
       su fail if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	



*Table with results
  #delimit ;
  esttab using "C:\Users\Fabio\Dropbox\py_mar_3\Analysis\output\reg_dur.tex", label replace 
  b(4) se(4) star(* 0.10 ** 0.05 *** 0.01) stats(mD N, labels("Dependent variable mean" "Observations") fmt(%9.0f) ) keep(tau) varlabels(tau Unilateral Divorce}) eqlabels(None) nonotes depvars nomtitles nonum
  
  prehead(\begin{table}[H]\centering 
				\scriptsize
				\caption{The average effect of unilateral divorce on the probability that a cohabitation spell ends (because of marriage or breakup). Unit of observation: cohabitation-month}  
				\label{tab:tabdur}
				\resizebox{0.8\textwidth}{!}{
				\begin{tabular}{l*{@M}{c}} \toprule
				\textit{Dependent variable:}&\multicolumn{4}{c}{\textit{keep cohabiting (0), cohabitation ends (1)}}\\)	 
			 posthead(\midrule)
             prefoot("State FE     & Yes & Yes& Yes & Yes\\
                      Year cohabitation began FE     & Yes & Yes& Yes & Yes\\
	                  Additional controls    & No & Yes& No& Yes\\
					  State-specific linear time trend    & No & No & Yes & Yes\\")   
			postfoot(
				\bottomrule \noalign{\smallskip}
				\end{tabular}
				}
				\begin{minipage}{\textwidth}
				\scriptsize\smallskip "Notes: This table reports the average effect of unilateral divorce on the probability of cohabitation ending in a breakup or marriage. The analysis follows the methodology outlined in \cite{borusyak2021} uses and the \textit{cohabitations sample}. The unit of observation is one cohabitation-month, which corresponds to a specific month within a particular cohabitation spell. The focus of this table is to report the probability of a cohabitation spell (initiated in state \textit{s} and year \textit{t}) ending in a breakup or marriage. The dummy variable \textit{Unilateral Divorce} takes value 1 if unilateral divorce was in effect in state \textit{s} and year \textit{t} at the start of the cohabitation spell. The additional controls include dummies for ethnicity, age, education, duration of cohabitation, and the introduction of unilateral divorce \textit{after} the cohabitation period commenced. Standard errors are clustered at the state level. Coefficients that are significantly different from zero are denoted by *10\%, **5\%, and ***1\%."  \\ 
				\end{minipage}
				\end{table}
				)
				  ;
	#delimit cr
	


*Event studies
did_imputation failc id  date_rel unil2 [weight=SAMWT] ,            fe($controls_weak ) autosample cluster(state1) maxit(200) allhorizons pre(8)
estimates store dur1
di e(pre_p)
did_imputation failc id  date_rel unil2 [weight=SAMWT] ,            fe($controls )      autosample cluster(state1) maxit(200) allhorizons pre(8)
estimates store dur2
di e(pre_p)
did_imputation failc id  date_rel unil2 [weight=SAMWT]  ,           fe($controls_weak ) autosample cluster(state1) maxit(200) allhorizons pre(8) timecontrols(state1) 
estimates store dur3
di e(pre_p)
did_imputation failc id  date_rel unil2 [weight=SAMWT]  ,           fe($controls )      autosample cluster(state1) maxit(200) allhorizons pre(8) timecontrols(state1) 
estimates store dur4
di e(pre_p)


#delimit ;
event_plot dur1 dur2 dur3 dur4, 
     plottype(scatter) ciplottype(rcap) 
	together perturb(-0.325(0.13)0.325) trimlead(7) trimlag(8) noautolegend 
	graph_opt(
		xtitle("Years since the event") ytitle("Average causal effect")  yscale(titlegap(*-25))
		xlabel(-7 "(-14,-13)" -6 "(-12,-11)" -5 "(-10,-9)" -4 "(-8,-7)" -3 "(-6,-5)" -2 "(-4,-3)" -1 "(-2,-1)"
		        0 "(0,1)"      1 "(2,3)"     2 "(4,5)"    3 "(6,7)"    4 "(8,9)"    5 "(10,11)"  6 "(10,11)"  7 "(12,13)" 8 "(14,15)" , format(%18.0fc) labsize(6pt)) ylabel(-0.04(0.01)0.04) 
		legend(order(1 "Basic controls" 3 "Full controls" 5 "Basic controls  + State-specific linear time trends" 7 "Full controls + State-specific linear time trends") rows(4) region(style(none))) 
		xline(-0.5, lcolor(gs8) lpattern(dash)) yline(0, lcolor(gs8)) graphregion(color(white)) bgcolor(white) ylabel(, angle(horizontal)) 
	) 
	lag_opt1(msymbol(Oh) color(cranberry)) lag_ci_opt1(color(cranberry)) 
	lag_opt2(msymbol(Sh) color(forest_green)) lag_ci_opt2(color(forest_green)) 	
	lag_opt3(msymbol(Dh) color(navy)) lag_ci_opt3(color(navy))
	lag_opt4(msymbol(Th) color(dkorange)) lag_ci_opt4(color(dkorange));
#delimit cr
graph export "C:\Users\Fabio\Dropbox\py_mar_3\Analysis\output\event_graph_dur.png",replace

/*
forvalues i=1(1)4{
#delimit ;
event_plot dur`i', default_look graph_opt(xtitle("Years since the event") ytitle("Coefficients") legend(off) 
 xlabel(-7 "(-14,-13)" -6 "(-12,-11)" -5 "(-10,-9)" -4 "(-8,-7)" -3 "(-6,-5)" -2 "(-4,-3)" -1 "(-2,-1)"
         0 "(0,1)"      1 "(2,3)"     2 "(4,5)"    3 "(6,7)"    4 "(8,9)"    5 "(10,11)"  6 "(10,11)"  7 "(12,13)" 8 "(14,15)" , 
		 format(%13.0fc) labsize(4pt)))  trimlead(7) trimlag(8);
#delimit cr
graph export "C:\Users\Fabio\Dropbox\py_mar_3\Analysis\output\event_graph_dur`i'.png",replace
}
*/
*******************************************
**3 Cohabitation Duration: risk of breakup
*******************************************
*Regressions
eststo clear
eststo:did_imputation failc id  date_rel_year unil2_year [weight=SAMWT],fe($controls_weak_year ) autosample cluster(state1)  
       su failc if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
	   
eststo:did_imputation failc id  date_rel_year unil2_year [weight=SAMWT],fe($controls_year )      autosample cluster(state1) 
       su failc if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	

eststo:did_imputation failc id  date_rel_year unil2_year [weight=SAMWT],fe($controls_weak_year ) autosample cluster(state1) timecontrols(state1) 
       su failc if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
	
eststo:did_imputation failc id  date_rel_year unil2_year [weight=SAMWT],fe($controls_year )      autosample cluster(state1) timecontrols(state1)   
       su failc if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	



*Table with results
  #delimit ;
  esttab using "C:\Users\Fabio\Dropbox\py_mar_3\Analysis\output\reg_dur_b.tex", label replace 
  b(4) se(4) star(* 0.10 ** 0.05 *** 0.01) stats(mD N, labels("Dependent variable mean" "Observations") fmt(%9.0f) ) keep(tau) varlabels(tau Unilateral Divorce}) eqlabels(None) nonotes depvars nomtitles nonum
  
  prehead(\begin{table}[H]\centering 
				\scriptsize
				\caption{The average effect of unilateral divorce on the probability that a cohabitation spell ends in a breakup. Unit of observation: cohabitation-month}  
				\label{tab:tabdurb}
				\resizebox{0.8\textwidth}{!}{
				\begin{tabular}{l*{@M}{c}} \toprule
				&\multicolumn{4}{c}{\textit{Dependent variable: keep cohabiting (0), breakup (1)}}\\)	 
			 posthead(\midrule)
             prefoot("State FE     & Yes & Yes& Yes & Yes\\
                      Year cohabitation began FE     & Yes & Yes& Yes & Yes\\
	                  Additional controls    & No & Yes& No& Yes\\
					  State-specific linear time trend    & No & No & Yes & Yes\\")   
			postfoot(
				\bottomrule \noalign{\smallskip}
				\end{tabular}
				}
				\begin{minipage}{\textwidth}
				\scriptsize\smallskip  "Notes: This table reports the average effect of unilateral divorce on the probability of cohabitation ending in a breakup. The analysis follows the methodology outlined in \cite{borusyak2021} uses and the \textit{cohabitations sample}. The unit of observation is one cohabitation-month, which corresponds to a specific month within a particular cohabitation spell. The focus of this table is to report the probability of a cohabitation spell (initiated in state \textit{s} and year \textit{t}) ending in a breakup, with the occurrence of marriage considered as a termination point (right censoring) for the cohabitation period. The dummy variable \textit{Unilateral Divorce} takes value 1 if unilateral divorce was in effect in state \textit{s} and year \textit{t} at the start of the cohabitation spell. The additional controls include dummies for ethnicity, age, education duration of cohabitation, and the introduction of unilateral divorce \textit{after} the cohabitation period commenced. Standard errors are clustered at the state level. Coefficients that are significantly different from zero are denoted by *10\%, **5\%, and ***1\%."  \\ 
				\end{minipage}\vspace{-6mm}
				\end{table}
				)
				  ;
	#delimit cr
	
	
*******************************************
**4 Cohabitation Duration: risk of marrying
*******************************************
*Regressions
eststo clear
eststo:did_imputation failm id  date_rel_year unil2_year [weight=SAMWT],fe($controls_weak_year ) autosample cluster(state1)  
       su failm if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
	   
eststo:did_imputation failm id  date_rel_year unil2_year [weight=SAMWT],fe($controls_year )      autosample cluster(state1) 
       su failm if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	

eststo:did_imputation failm id  date_rel_year unil2_year [weight=SAMWT],fe($controls_weak_year ) autosample cluster(state1) timecontrols(state1) 
       su failm if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
	
eststo:did_imputation failm id  date_rel_year unil2_year [weight=SAMWT],fe($controls_year )      autosample cluster(state1) timecontrols(state1)   
       su failm if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	



*Table with results
  #delimit ;
  esttab using "C:\Users\Fabio\Dropbox\py_mar_3\Analysis\output\reg_dur_m.tex", label replace 
  b(4) se(4) star(* 0.10 ** 0.05 *** 0.01) stats(mD N, labels("Dependent variable mean" "Observations") fmt(%9.0f) ) keep(tau) varlabels(tau Unilateral Divorce}) eqlabels(None) nonotes depvars nomtitles nonum
  
  prehead(\begin{table}[H]\centering 
				\scriptsize
				\caption{The average effect of unilateral divorce on the probability that a cohabitation spell ends in a marriage. Unit of observation: cohabitation-month}  
				\label{tab:tabdurm}
				\resizebox{0.8\textwidth}{!}{
				\begin{tabular}{l*{@M}{c}} \toprule
				&\multicolumn{4}{c}{\textit{Dependent variable: keep cohabiting (0), marry (1)}}\\)	 
			 posthead(\midrule)
             prefoot("State FE     & Yes & Yes& Yes & Yes\\
                      Year cohabitation began FE     & Yes & Yes& Yes & Yes\\
	                  Additional controls    & No & Yes& No& Yes\\
					  State-specific linear time trend    & No & No & Yes & Yes\\")   
			postfoot(
				\bottomrule \noalign{\smallskip}
				\end{tabular}
				}
				\begin{minipage}{\textwidth}
				\scriptsize\smallskip "Notes: This table reports the average effect of unilateral divorce on the probability of cohabitation ending in a marriage. The analysis follows the methodology outlined in \cite{borusyak2021} uses and the \textit{cohabitations sample}. The unit of observation is one cohabitation-month, which corresponds to a specific month within a particular cohabitation spell. The focus of this table is to report the probability of a cohabitation spell (initiated in state \textit{s} and year \textit{t}) ending in a marriage, with the occurrence of breakup considered as a termination point (right censoring) for the cohabitation period. The dummy variable \textit{Unilateral Divorce} takes value 1 if unilateral divorce was in effect in state \textit{s} and year \textit{t} at the start of the cohabitation spell. The additional controls include dummies for ethnicity, age, education, duration of cohabitation, and the introduction of unilateral divorce \textit{after} the cohabitation period commenced. Standard errors are clustered at the state level. Coefficients that are significantly different from zero are denoted by *10\%, **5\%, and ***1\%."  \\ 
				\end{minipage}
				\end{table}
				)
				  ;
	#delimit cr
	


	
******************************************
**5 Robustness: Only community property breakup
*******************************************

*Regressions
eststo clear
eststo:did_imputation failc id  date_rel_year unil2_year [weight=SAMWT] if com==1,fe($controls_weak_year ) autosample cluster(state1)  maxit(400)
       su fail if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
	   
eststo:did_imputation failc id  date_rel_year unil2_year [weight=SAMWT]  if com==1,fe($controls_year )      autosample cluster(state1) maxit(400)
       su fail if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	

eststo:did_imputation failc id  date_rel_year unil2_year [weight=SAMWT] if tit==1,fe($controls_weak_year ) autosample cluster(state1)  maxit(400)
       su fail if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
	   
eststo:did_imputation failc id  date_rel_year unil2_year [weight=SAMWT]  if tit==1,fe($controls_year )      autosample cluster(state1) maxit(400)
       su fail if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	



*Table with results
  #delimit ;
  esttab using "C:\Users\Fabio\Dropbox\py_mar_3\Analysis\output\reg_dur_b_com.tex", label replace 
  b(4) se(4) star(* 0.10 ** 0.05 *** 0.01) stats(mD N, labels("Dependent variable mean" "Observations") fmt(%9.0f) ) keep(tau) varlabels(tau Unilateral Divorce}) eqlabels(None) nonotes depvars nomtitles nonum
  
  prehead(\begin{table}[H]\centering 
				\scriptsize
				\caption{The average effect of unilateral divorce on the probability that a cohabitation spell ends in a breakup by property regime at divorce. Unit of observation: cohabitation-month}  
				\label{tab:tabdurbcom}
				\resizebox{0.8\textwidth}{!}{
				\begin{tabular}{l*{@M}{c}} \toprule
				&\multicolumn{4}{c}{\textit{Dependent variable: keep cohabiting (0), breakup (1)}}\\
				\textit{Sample:}& ComP & ComP &  Tit & Tit \\)		 			
			 posthead(\midrule)
             prefoot("State FE     & Yes & Yes& Yes & Yes\\
                      Year cohabitation began FE     & Yes & Yes& Yes & Yes\\
	                  Additional controls    & No & Yes& No& Yes\\")   
			postfoot(
				\bottomrule \noalign{\smallskip}
				\end{tabular}
				}
				\begin{minipage}{\textwidth}
				\scriptsize\smallskip  "Notes: This table reports the average effect of unilateral divorce on the probability of cohabitation ending in a breakup. The analysis follows the methodology outlined in \cite{borusyak2021} uses and the \textit{cohabitations sample}. The unit of observation is one cohabitation-month, which corresponds to a specific month within a particular cohabitation spell. The focus of this table is to report the probability of a cohabitation spell (initiated in state \textit{s} and year \textit{t}) ending in a breakup, with the occurrence of marriage considered as a termination point (right censoring) for the cohabitation period. The dummy variable \textit{Unilateral Divorce} takes value 1 if unilateral divorce was in effect in state \textit{s} and year \textit{t} at the start of the cohabitation spell. The additional controls include dummies for ethnicity, age, education, duration of cohabitation,  and the introduction of unilateral divorce \textit{after} the cohabitation period commenced. ComP: community property states; Tit: title based regime. Standard errors are clustered at the state level. Coefficients that are significantly different from zero are denoted by *10\%, **5\%, and ***1\%."  \\ 
				\end{minipage}
				\end{table}
				)
				  ;
	#delimit cr
	
	
	
******************************************
**6 Robustness: Only community property marriage
*******************************************

*Regressions
eststo clear
eststo:did_imputation failm id  date_rel_year unil2_year [weight=SAMWT] if com==1,fe($controls_weak_year ) autosample cluster(state1)  maxit(400)
       su fail if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
	   
eststo:did_imputation failm id  date_rel_year unil2_year [weight=SAMWT]  if com==1,fe($controls_year )      autosample cluster(state1) maxit(400)
       su fail if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	

eststo:did_imputation failm id  date_rel_year unil2_year [weight=SAMWT] if tit==1,fe($controls_weak_year ) autosample cluster(state1)  maxit(400)
       su fail if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
	   
eststo:did_imputation failm id  date_rel_year unil2_year [weight=SAMWT]  if tit==1,fe($controls_year )      autosample cluster(state1) maxit(400)
       su fail if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	



*Table with results
  #delimit ;
  esttab using "C:\Users\Fabio\Dropbox\py_mar_3\Analysis\output\reg_dur_m_com.tex", label replace 
  b(4) se(4) star(* 0.10 ** 0.05 *** 0.01) stats(mD N, labels("Dependent variable mean" "Observations") fmt(%9.0f) ) keep(tau) varlabels(tau Unilateral Divorce}) eqlabels(None) nonotes depvars nomtitles nonum
  
  prehead(\begin{table}[H]\centering 
				\scriptsize
				\caption{The average effect of unilateral divorce on the probability that a cohabitation spell ends in a marriage by property regime at divorce. Unit of observation: cohabitation-month}  
				\label{tab:tabdurmcom}
				\resizebox{0.8\textwidth}{!}{
				\begin{tabular}{l*{@M}{c}} \toprule
				&\multicolumn{4}{c}{\textit{Dependent variable: keep cohabiting (0), marry (1)}}\\
				\textit{Sample:}& ComP & ComP &  Tit & Tit \\)		 
			 posthead(\midrule) 
             prefoot("State FE     & Yes & Yes& Yes & Yes\\
                      Year cohabitation began FE     & Yes & Yes& Yes & Yes\\
	                  Additional controls    & No & Yes& No& Yes\\")   
			postfoot(
				\bottomrule \noalign{\smallskip}
				\end{tabular}
				}
				\begin{minipage}{\textwidth}
				\scriptsize\smallskip  "Notes: This table reports the average effect of unilateral divorce on the probability of cohabitation ending in a marriage. The analysis follows the methodology outlined in \cite{borusyak2021} uses and the \textit{cohabitations sample}. The unit of observation is one cohabitation-month, which corresponds to a specific month within a particular cohabitation spell. The focus of this table is to report the probability of a cohabitation spell (initiated in state \textit{s} and year \textit{t}) ending in a marriage, with the occurrence of breakup considered as a termination point (right censoring) for the cohabitation period. The dummy variable \textit{Unilateral Divorce} takes value 1 if unilateral divorce was in effect in state \textit{s} and year \textit{t} at the start of the cohabitation spell. The additional controls include dummies for ethnicity, age, education, duration of cohabitation, and the introduction of unilateral divorce \textit{after} the cohabitation period commenced. ComP: community property states; Tit: title based regime. Standard errors are clustered at the state level. Coefficients that are significantly different from zero are denoted by *10\%, **5\%, and ***1\%."  \\ 
				\end{minipage}
				\end{table}
				)
				  ;
	#delimit cr
	

*******************************************
**7 Cohabitation Duration: risk of breakup only residents
*******************************************
*Regressions
eststo clear
eststo:did_imputation failc id  date_rel_year unil2_year [weight=SAMWT] if keep==1,fe($controls_weak_year ) autosample cluster(state1)  
       su failc if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
	   
eststo:did_imputation failc id  date_rel_year unil2_year [weight=SAMWT] if keep==1,fe($controls_year )      autosample cluster(state1) 
       su failc if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	

eststo:did_imputation failc id  date_rel_year unil2_year [weight=SAMWT] if keep==1,fe($controls_weak_year ) autosample cluster(state1) timecontrols(state1) 
       su failc if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
	
eststo:did_imputation failc id  date_rel_year unil2_year [weight=SAMWT] if keep==1,fe($controls_year )      autosample cluster(state1) timecontrols(state1)   
       su failc if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	



*Table with results
  #delimit ;
  esttab using "C:\Users\Fabio\Dropbox\py_mar_3\Analysis\output\reg_dur_b_keep.tex", label replace 
  b(4) se(4) star(* 0.10 ** 0.05 *** 0.01) stats(mD N, labels("Dependent variable mean" "Observations") fmt(%9.0f) ) keep(tau) varlabels(tau Unilateral Divorce}) eqlabels(None) nonotes depvars nomtitles nonum
  
  prehead(\begin{table}[H]\centering 
				\scriptsize
				\caption{The average effect of unilateral divorce on the probability that a cohabitation spell ends in a breakup. Unit of observation: cohabitation-month. Sample of never movers}  
				\label{tab:tabdurbk}
				\resizebox{0.8\textwidth}{!}{
				\begin{tabular}{l*{@M}{c}} \toprule
				&\multicolumn{4}{c}{\textit{Dependent variable: keep cohabiting (0), breakup (1)}}\\)	 
			 posthead(\midrule)
             prefoot("State FE     & Yes & Yes& Yes & Yes\\
                      Year cohabitation began FE     & Yes & Yes& Yes & Yes\\
	                  Additional controls    & No & Yes& No& Yes\\
					  State-specific linear time trend    & No & No & Yes & Yes\\")   
			postfoot(
				\bottomrule \noalign{\smallskip}
				\end{tabular}
				}
				\begin{minipage}{\textwidth}
				\scriptsize\smallskip  "Notes: This table reports the average effect of unilateral divorce on the probability of cohabitation ending in a breakup. The analysis follows the methodology outlined in \cite{borusyak2021} uses and the \textit{cohabitations sample}. The unit of observation is one cohabitation-month, which corresponds to a specific month within a particular cohabitation spell.  The sample used for this regression includes only respondents who still live in the same state where they were living at sixteen years old. The focus of this table is to report the probability of a cohabitation spell (initiated in state \textit{s} and year \textit{t}) ending in a breakup, with the occurrence of marriage considered as a termination point (right censoring) for the cohabitation period. The dummy variable \textit{Unilateral Divorce} takes value 1 if unilateral divorce was in effect in state \textit{s} and year \textit{t} at the start of the cohabitation spell. The additional controls include dummies for ethnicity, age, education, duration of cohabitation, and the introduction of unilateral divorce \textit{after} the cohabitation period commenced. Standard errors are clustered at the state level. Coefficients that are significantly different from zero are denoted by *10\%, **5\%, and ***1\%."  \\ 
				\end{minipage}
				\end{table}
				)
				  ;
	#delimit cr
	
	
*******************************************
**8 Cohabitation Duration: risk of marrying only residents
*******************************************
*Regressions
eststo clear
eststo:did_imputation failm id  date_rel_year unil2_year [weight=SAMWT] if keep==1,fe($controls_weak_year ) autosample cluster(state1)  
       su failm if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
	   
eststo:did_imputation failm id  date_rel_year unil2_year [weight=SAMWT] if keep==1,fe($controls_year )      autosample cluster(state1) 
       su failm if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	

eststo:did_imputation failm id  date_rel_year unil2_year [weight=SAMWT] if keep==1,fe($controls_weak_year ) autosample cluster(state1) timecontrols(state1) 
       su failm if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
	
eststo:did_imputation failm id  date_rel_year unil2_year [weight=SAMWT] if keep==1,fe($controls_year )      autosample cluster(state1) timecontrols(state1)   
       su failm if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	



*Table with results
  #delimit ;
  esttab using "C:\Users\Fabio\Dropbox\py_mar_3\Analysis\output\reg_dur_m_keep.tex", label replace 
  b(4) se(4) star(* 0.10 ** 0.05 *** 0.01) stats(mD N, labels("Dependent variable mean" "Observations") fmt(%9.0f) ) keep(tau) varlabels(tau Unilateral Divorce}) eqlabels(None) nonotes depvars nomtitles nonum
  
  prehead(\begin{table}[H]\centering 
				\scriptsize
				\caption{The average effect of unilateral divorce on the probability that a cohabitation spell ends in a marriage. Unit of observation: cohabitation-month. Sample of never movers}  
				\label{tab:tabdurmk}
				\resizebox{0.8\textwidth}{!}{
				\begin{tabular}{l*{@M}{c}} \toprule
				&\multicolumn{4}{c}{\textit{Dependent variable: keep cohabiting (0), marry (1)}}\\)	 
			 posthead(\midrule)
             prefoot("State FE     & Yes & Yes& Yes & Yes\\
                      Year cohabitation began FE     & Yes & Yes& Yes & Yes\\
	                  Additional controls    & No & Yes& No& Yes\\
					  State-specific linear time trend    & No & No & Yes & Yes\\")   
			postfoot(
				\bottomrule \noalign{\smallskip}
				\end{tabular}
				}
				\begin{minipage}{\textwidth}
				\scriptsize\smallskip "Notes: This table reports the average effect of unilateral divorce on the probability of cohabitation ending in a marriage. The analysis follows the methodology outlined in \cite{borusyak2021} uses and the \textit{cohabitations sample}. The unit of observation is one cohabitation-month, which corresponds to a specific month within a particular cohabitation spell.  The sample used for this regression includes only respondents who still live in the same state where they were living at sixteen years old. The focus of this table is to report the probability of a cohabitation spell (initiated in state \textit{s} and year \textit{t}) ending in a marriage, with the occurrence of breakup considered as a termination point (right censoring) for the cohabitation period. The dummy variable \textit{Unilateral Divorce} takes value 1 if unilateral divorce was in effect in state \textit{s} and year \textit{t} at the start of the cohabitation spell. The additional controls include dummies for ethnicity, age, education, duration of cohabitation, and the introduction of unilateral divorce \textit{after} the cohabitation period commenced. Standard errors are clustered at the state level. Coefficients that are significantly different from zero are denoted by *10\%, **5\%, and ***1\%."  \\ 
				\end{minipage}
				\end{table}
				)
				  ;
	#delimit cr