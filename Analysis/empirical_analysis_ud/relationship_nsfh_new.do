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

****************************************
*BUILD SECOND RELATIONSHIP HERE
***************************************
gen inspect=0
replace inspect=1 if NUMCOHAB+M95>=2

*First Marriage End
gen broken_m=0 if mar==1
replace broken_m=1 if mar==1 & M99!=6

*First Cohabitation End
gen brocken_c=0 if coh==1
*Separation
replace brocken_c=1 if coh==1 & HOWEND1==1
*Marriage
replace brocken_c=2 if coh==1 & HOWEND1>=2 & HOWEND1<=3  
*No Information-dont check this
replace inspect=0 if HOWEND1==9 & coh==1

*Dont consider cases where 1 m+1c and transition from one to the other
replace inspect=0 if NUMCOHAB+M95==2 & brocken_c==2

*Now create second relationship
gen dou=1
replace dou=2 if inspect==1

*Expand the data
expand dou, generate(reln)
replace reln=reln+1

*Get info on second relationship

*Get second marriage if already married
replace date_m=M103T02M if mar==1 & reln==2 & M103T02M!=9996
replace date_c=BEGCOH2 if coh==1 & reln==2 & BEGCOH2!=0

replace mar=. if reln==2
replace coh=. if reln==2

replace mar=1 if date_m!=. & date_c>=. & reln==2
replace coh=0 if date_m!=. & date_c>=. & reln==2

replace mar=0 if date_m>=. & date_c!=. & reln==2
replace coh=1 if date_m>=. & date_c!=. & reln==2

replace mar=1 if date_m<=date_c & date_m!=. & date_c!=. & reln==2
replace coh=0 if date_m<=date_c & date_m!=. & date_c!=. & reln==2

replace mar=0 if date_m>date_c & date_m!=. & date_c!=. & reln==2
replace coh=1 if date_m>date_c & date_m!=. & date_c!=. & reln==2

replace date_rel=date_m if mar==1 & reln==2
replace date_rel=date_c if mar==0 & reln==2
drop if date_rel>2000 & reln==2

replace age_rel=(date_rel-birth)/12 if reln==2

*replace clength=BEGCOH2-ENDCOH2 if clength==. &

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
drop if clength<=0 & mar==0 

*Keep if not too far away from policy change
*keep if (((date_rel-unill<=20) & (date_rel-unill>=-20)) | unill==0 )
*keep if (((date_rel-unill>=-20)) | unill==0 )


*What to do with age at first relationship<16?
*Many people fall here, check!
keep if date_rel>=1950 


********************************************************************
*E) Controls
********************************************************************

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

*Time for unilateral divorce
gen unil2=unill
replace unil2=. if unill==0

*Etnicity
gen etn=0
replace etn=1 if M484==1
replace etn=2 if M484==2
replace etn=3 if M484>=3

*Month at which relationship started 
gen monthh=date_relo-(date_rel-1900)*12

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

gen agem=floor(M2BP01/10)

*Save file for doing fancy graphs with matplotlib
gen years=floor((date_rel_year-unil2_year)/1)*1
gen date_rel_agg=floor(date_rel_year/1)*1
gen wgt=SAMWT

savesome mar unidd years date_rel_agg wgt using "C:\Users\Fabio\Dropbox\py_mar_3\Analysis\output\fig1.dta", replace
********************************************************************************
********************************************************************************
***Actual analysis of divorce laws
********************************************************************************
********************************************************************************

**********************
**0 Summary statistics
*********************
gen cohabitt=0
replace cohabitt=1 if mar==0


label variable mar "Marry (0/1)"
label variable cohabitt "Cohabit (0/1)"
label variable M2BP01 "Age of respondent"
label variable coll "College (0/1)"
label variable date_rel_year "Year the Cohabitation began"
label variable unidd "Unilateral divorce when relationship began (0/1)"


est clear
estpost tabstat mar cohabitt  coll date_rel_year unidd M2BP01 [weight=SAMWT], c(stat) stat(count mean median sd)
#delimit ;
esttab using "C:\Users\Fabio\Dropbox\py_mar_3\Analysis\output\sum_sing.tex", replace
 cells("count mean(fmt(%6.2fc)) p50 sd(fmt(%6.2fc))")   nonumber 
  nomtitle nonote noobs label booktabs 
  collabels("N" "Mean" "Median" "SD") ;
#delimit cr


********************************************************************************
*1) Basic regression
********************************************************************************

*Create controls
global controls      monthh state1 date_rel coll etn agem
global controls_weak state1 date_rel 
global controls_year   monthh state1 date_rel_year coll etn agem
global controls_weak_year state1 date_rel_year


*Regressions
eststo clear
eststo:did_imputation mar id  date_rel_year unil2_year [weight=SAMWT],fe($controls_weak_year ) autosample cluster(state1) 
       su mar [weight=SAMWT] if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
	   
eststo:did_imputation mar id  date_rel_year unil2_year [weight=SAMWT],fe($controls_year )      autosample cluster(state1)  
       su mar [weight=SAMWT]  if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
eststo:did_imputation mar id  date_rel_year unil2_year [weight=SAMWT],fe($controls_weak_year ) autosample cluster(state1)  timecontrols(state1) 
       su mar [weight=SAMWT]  if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
eststo:did_imputation mar id  date_rel_year unil2_year [weight=SAMWT],fe($controls_year )      autosample cluster(state1)  timecontrols(state1) 
       su mar [weight=SAMWT]  if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
	   


*Table with results
  #delimit ;
  esttab using "C:\Users\Fabio\Dropbox\py_mar_3\Analysis\output\reg_rel.tex", label replace 
    b(4) se(4) star(* 0.10 ** 0.05 *** 0.01) stats(mD N, labels("Dependent variable mean" "Observations") fmt(%9.0f) ) keep(tau) varlabels(tau Unilateral Divorce}) eqlabels(None) nonotes depvars nomtitles nonum
  
  prehead(\begin{table}[H]\centering 
				\scriptsize
				\caption{The average effect of unilateral divorce on the choice between cohabitation and marriage among newly formed couples}  
				\label{tab:tabrel}
				\resizebox{0.8\textwidth}{!}{
				\begin{tabular}{l*{@M}{c}} \toprule
				&\multicolumn{4}{c}{\textit{Dependent variable: marry (1), cohabit (0)} }\\)	 
			 posthead(\midrule)
             prefoot("State FE     & Yes & Yes& Yes & Yes\\
                      Period relationship began FE     & Yes & Yes& Yes & Yes\\
	                  Additional controls    & No & Yes& No& Yes\\
					  State-specific linear time trend    & No & No & Yes & Yes\\")   
			postfoot(
				\bottomrule \noalign{\smallskip}
				\end{tabular}
				}
				\begin{minipage}{\textwidth}
				\scriptsize\smallskip "Notes: This table reports the average effect of unilateral divorce on the choice between cohabitation and marriage. The analysis follows the methodology outlined in \cite{borusyak2021} uses and the \textit{relationships sample}. The unit of observation refers to newly formed couples (either married or cohabiting) in state \textit{s} and year \textit{t}. The dummy variable \textit{Unilateral Divorce} takes the value 1 if unilateral divorce was in effect in state \textit{s} and year \textit{t} and 0 otherwise. The additional controls include dummies for ethnicity, age, and education. Standard errors are clustered at the state level. Coefficients that are significantly different from zero are denoted by *10\%, **5\%, and ***1\%."   \\ 
				\end{minipage}\vspace{-6mm}
				\end{table}
				)
				  ;
	#delimit cr


*Event studies
did_imputation mar id  date_rel unil2 [weight=SAMWT] ,            fe($controls_weak ) autosample cluster(state1) maxit(200) allhorizons pre(8)
estimates store dur1
di e(pre_p)
did_imputation mar id  date_rel unil2 [weight=SAMWT] ,            fe($controls )      autosample cluster(state1) maxit(200) allhorizons pre(8)
estimates store dur2
di e(pre_p)
did_imputation mar id  date_rel unil2 [weight=SAMWT]  ,           fe($controls_weak ) autosample cluster(state1) maxit(200) allhorizons pre(8) timecontrols(state1) 
estimates store dur3
di e(pre_p)
did_imputation mar id  date_rel unil2 [weight=SAMWT]  ,           fe($controls )      autosample cluster(state1) maxit(200) allhorizons pre(8) timecontrols(state1) 
estimates store dur4
di e(pre_p)


#delimit ;
event_plot dur1 dur2 dur3 dur4, 
     plottype(scatter) ciplottype(rcap) 
	together perturb(-0.325(0.13)0.325) trimlead(7) trimlag(8) noautolegend 
	graph_opt(
		xtitle("Years since the event") ytitle("Average causal effect")  yscale(titlegap(*-30))
		xlabel(-7 "(-14,-13)" -6 "(-12,-11)" -5 "(-10,-9)" -4 "(-8,-7)" -3 "(-6,-5)" -2 "(-4,-3)" -1 "(-2,-1)"
		        0 "(0,1)"      1 "(2,3)"     2 "(4,5)"    3 "(6,7)"    4 "(8,9)"    5 "(10,11)"  6 "(10,11)"  7 "(12,13)" 8 "(14,15)" , format(%18.0fc) labsize(6pt)) ylabel(-0.3(0.1)0.3) 
		legend(order(1 "Basic controls" 3 "Full controls" 5 "Basic controls  + State-specific linear time trends" 7 "Full controls + State-specific linear time trends") rows(4) region(style(none))) 
		xline(-0.5, lcolor(gs8) lpattern(dash)) yline(0, lcolor(gs8)) graphregion(color(white)) bgcolor(white) ylabel(, angle(horizontal)) 
	) 
	lag_opt1(msymbol(Oh) color(cranberry)) lag_ci_opt1(color(cranberry)) 
	lag_opt2(msymbol(Sh) color(forest_green)) lag_ci_opt2(color(forest_green)) 	
	lag_opt3(msymbol(Dh) color(navy)) lag_ci_opt3(color(navy))
	lag_opt4(msymbol(Th) color(dkorange)) lag_ci_opt4(color(dkorange));
#delimit cr
graph export "C:\Users\Fabio\Dropbox\py_mar_3\Analysis\output\event_graph_rel.png",replace

/*
forvalues i=1(1)4{
#delimit ;
event_plot dur`i', default_look graph_opt(xtitle("Years since the event") ytitle("Coefficients") legend(off) 
 xlabel(-7 "(-14,-13)" -6 "(-12,-11)" -5 "(-10,-9)" -4 "(-8,-7)" -3 "(-6,-5)" -2 "(-4,-3)" -1 "(-2,-1)"
         0 "(0,1)"      1 "(2,3)"     2 "(4,5)"    3 "(6,7)"    4 "(8,9)"    5 "(10,11)"  6 "(10,11)"  7 "(12,13)" 8 "(14,15)" , 
		 format(%18.0fc) labsize(4pt)))  trimlead(7) trimlag(8);
#delimit cr
graph export "C:\Users\Fabio\Dropbox\py_mar_3\Analysis\output\event_graph_rel`i'.png",replace
}
*/


********************************************************************************
*2) By Divorce property regime
********************************************************************************

eststo clear
eststo:did_imputation mar id  date_rel_year unil2_year [weight=SAMWT] if com==1,fe($controls_weak_year )      autosample cluster(state1)  
       su mar [weight=SAMWT]  if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
	   scalar beta_c1= e(b)[1,1]
	   scalar v_beta_c1= e(V)[1,1]
	   



eststo:did_imputation mar id  date_rel_year unil2_year [weight=SAMWT] if com==1,fe($controls_year )      autosample cluster(state1)  
       su mar [weight=SAMWT]  if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
	   scalar beta_c2= e(b)[1,1]
	   scalar v_beta_c2= e(V)[1,1]
	   
	   
eststo:did_imputation mar id  date_rel_year unil2_year [weight=SAMWT] if tit==1,fe($controls_weak_year )      autosample cluster(state1) 
       su mar [weight=SAMWT]  if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
	   scalar beta_t1= e(b)[1,1]
	   scalar v_beta_t1= e(V)[1,1]
	   
eststo:did_imputation mar id  date_rel_year unil2_year [weight=SAMWT] if tit==1,fe($controls_year )      autosample cluster(state1) 
       su mar [weight=SAMWT]  if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
	   scalar beta_t2= e(b)[1,1]
	   scalar v_beta_t2= e(V)[1,1]

	   
scalar pvalue1=normal((beta_c1-beta_t1)/(v_beta_c1+v_beta_t1)^0.5)
scalar pvalue2=normal((beta_c2-beta_t2)/(v_beta_c2+v_beta_t2)^0.5)
di pvalue1
di pvalue2




*Table with results
  #delimit ;
  esttab using "C:\Users\Fabio\Dropbox\py_mar_3\Analysis\output\reg_rel_com.tex", label replace 
  b(4) se(4) star(* 0.10 ** 0.05 *** 0.01) stats(mD N, labels("Dependent variable mean" "Observations") fmt(%9.0f) ) keep(tau) varlabels(tau Unilateral Divorce}) eqlabels(None) nonotes depvars nonum nomtitles
  
  prehead(\begin{table}[H]\centering 
				\scriptsize

				\caption{The average effect of unilateral divorce on the choice between cohabitation and marriage among newly formed couples, by property regime upon divorce}  
                 \label{tab:tabrelcom}
				\resizebox{0.8\textwidth}{!}{
				\begin{tabular}{l*{@M}{c}} \toprule
				&\multicolumn{4}{c}{\textit{Dependent variable: marry (1), cohabit (0)}}\\
				\textit{Sample:}& ComP & ComP &  Tit & Tit \\)	 
			 posthead(\midrule)
             prefoot("State FE     & Yes & Yes& Yes & Yes\\
                      Period relationship began FE     & Yes & Yes& Yes & Yes\\
	                  Additional controls    & No & Yes& No& Yes\\")   
			postfoot(
				\bottomrule \noalign{\smallskip}
				\end{tabular}
				}
				\begin{minipage}{\textwidth}
				\scriptsize\smallskip "Notes: This table reports the average effect of unilateral divorce on the choice between cohabitation and marriage. The analysis follows the methodology outlined in \cite{borusyak2021} uses and the \textit{relationships sample}. The unit of observation refers to newly formed couples (either married or cohabiting) in state \textit{s} and year \textit{t}. The dummy variable \textit{Unilateral Divorce} takes the value 1 if unilateral divorce was in effect in state \textit{s} and year \textit{t} and 0 otherwise. The additional controls include dummies for ethnicity, age, and education. ComP: community property states; Tit: title based regime. Standard errors are clustered at the state level. Coefficients that are significantly different from zero are denoted by *10\%, **5\%, and ***1\%."  \\ 
				\end{minipage}
				\end{table}
				)
				  ;
	#delimit cr

	
********************************************************************************
*3) Heterogeneity by Children
********************************************************************************

eststo clear
eststo:did_imputation mar id  date_rel_year unil2_year [weight=SAMWT] if nch>=1,fe($controls_year )      autosample cluster(state1)  
       su mar [weight=SAMWT]  if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	

	   


	   
	   
eststo:did_imputation mar id  date_rel_year unil2_year [weight=SAMWT] if nch==0 ,fe($controls_year )      autosample cluster(state1) 
       su mar [weight=SAMWT]  if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	

	   






*Table with results
  #delimit ;
  esttab using "C:\Users\Fabio\Dropbox\py_mar_3\Analysis\output\reg_rel_ch.tex", label replace 
  b(4) se(4) star(* 0.10 ** 0.05 *** 0.01) stats(mD N, labels("Dependent variable mean" "Observations") fmt(%9.0f) ) keep(tau) varlabels(tau Unilateral Divorce}) eqlabels(None) nonotes depvars nonum nomtitles
  
  prehead(\begin{table}[H]\centering 
				\scriptsize
				
				\caption{The average effect of unilateral divorce on the choice between cohabitation and marriage among newly formed couples, by the presence of children}  
				\label{tab:tabrelch}
				\resizebox{0.8\textwidth}{!}{
				\begin{tabular}{l*{@M}{c}} \toprule
				&\multicolumn{2}{c}{\textit{Dependent variable: marry (1), cohabit (0)}}\\
				\textit{Sample:}& Some children &  Childless \\)	 
			 posthead(\midrule)
             prefoot("State FE     & Yes & Yes \\
                      Period relationship began FE     & Yes & Yes \\
	                  Additional controls    & Yes & Yes\\")   
			postfoot(
				\bottomrule \noalign{\smallskip}
				\end{tabular}
				}
				\begin{minipage}{\textwidth}
				\scriptsize\smallskip "Notes: This table reports the average effect of unilateral divorce on the choice between cohabitation and marriage. The analysis follows the methodology outlined in \cite{borusyak2021} uses and the \textit{relationships sample}. The unit of observation refers to newly formed couples (either married or cohabiting) in state \textit{s} and year \textit{t}. The dummy variable \textit{Unilateral Divorce} takes the value 1 if unilateral divorce was in effect in state \textit{s} and year \textit{t} and 0 otherwise. The additional controls include dummies for ethnicity, age, and education. Standard errors are clustered at the state level. Coefficients that are significantly different from zero are denoted by *10\%, **5\%, and ***1\%."  \\ 
				\end{minipage}
				\end{table}
				)
				  ;
	#delimit cr


	
********************************************************************************
*4) Only Residents
********************************************************************************





*Regressions
eststo clear
eststo:did_imputation mar id  date_rel_year unil2_year [weight=SAMWT] if keep==1,fe($controls_weak_year ) autosample cluster(state1) 
       su mar [weight=SAMWT] if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
	   
eststo:did_imputation mar id  date_rel_year unil2_year [weight=SAMWT] if keep==1,fe($controls_year )      autosample cluster(state1)  
       su mar [weight=SAMWT]  if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
eststo:did_imputation mar id  date_rel_year unil2_year [weight=SAMWT] if keep==1,fe($controls_weak_year ) autosample cluster(state1)  timecontrols(state1) 
       su mar [weight=SAMWT]  if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
eststo:did_imputation mar id  date_rel_year unil2_year [weight=SAMWT] if keep==1,fe($controls_year )      autosample cluster(state1)  timecontrols(state1) 
       su mar [weight=SAMWT]  if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
	   


*Table with results
  #delimit ;
  esttab using "C:\Users\Fabio\Dropbox\py_mar_3\Analysis\output\reg_rel_keep.tex", label replace 
    b(4) se(4) star(* 0.10 ** 0.05 *** 0.01) stats(mD N, labels("Dependent variable mean" "Observations") fmt(%9.0f) ) keep(tau) varlabels(tau Unilateral Divorce}) eqlabels(None) nonotes depvars nomtitles nonum
  
  prehead(\begin{table}[H]\centering 
				\scriptsize
				
				\caption{The average effect of unilateral divorce on the choice of cohabitation vs. marriage among newly formed couples. Sample of never movers}  
				\label{tab:tabrelkeep}
				\resizebox{0.8\textwidth}{!}{
				\begin{tabular}{l*{@M}{c}} \toprule
				&\multicolumn{4}{c}{\textit{Dependent variable: marry (1), cohabit (0)}}\\)	 
			 posthead(\midrule)
             prefoot("State FE     & Yes & Yes& Yes & Yes\\
                      Period relationship began FE     & Yes & Yes& Yes & Yes\\
	                  Additional controls    & No & Yes& No& Yes\\
					  State-specific linear time trend    & No & No & Yes & Yes\\")   
			postfoot(
				\bottomrule \noalign{\smallskip}
				\end{tabular}
				}
				\begin{minipage}{\textwidth}
				\scriptsize\smallskip "Notes: This table reports the average effect of unilateral divorce on the choice between cohabitation and marriage. The analysis follows the methodology outlined in \cite{borusyak2021} uses and the \textit{relationships sample}. The sample used for this regression includes only respondents who still live in the same state where they were living at sixteen years old. The unit of observation refers to newly formed couples (either married or cohabiting) in state \textit{s} and year \textit{t}. The dummy variable \textit{Unilateral Divorce} takes the value 1 if unilateral divorce was in effect in state \textit{s} and year \textit{t} and 0 otherwise. The additional controls include dummies for ethnicity, age, and education. Standard errors are clustered at the state level. Coefficients that are significantly different from zero are denoted by *10\%, **5\%, and ***1\%."  \\ 
				\end{minipage}
				\end{table}
				)
				  ;
	#delimit cr

