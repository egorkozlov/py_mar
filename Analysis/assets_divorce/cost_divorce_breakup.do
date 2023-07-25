***************************
*Cost of Divorce in the PSID
*****************************



*********************************
*Get assets
*********************************
clear all
import delimited "C:\Users\Fabio\Dropbox\py_mar_3\Analysis\assets_divorce\assets.csv", asdouble 
quietly destring *, replace ignore("NA")



global clean house cash bonds stocks penval carval other_debt  mort1 mort2 real_estate busval



foreach val in $clean2{

replace `val'=. if `val'>=999999998
replace `val'=. if `val'<=-99999998
                          
}
replace mort1=. if mort1>=9999998
replace mort2=. if mort2>=9999998
replace stocks=. if stocks>=99999998

*Generate variations of wealth
g wealth=house+cash+bonds+stocks+penval+carval-other_debt -mort1-mort2+real_estate+busval

*Compare our constructed wealth measure with the PSID one for 2009: similar variable and values!
gen genera=a1+a2+a3+a4+a5+a6+a7-d1
corr genera constructed_net_worth if year==2009 & genera>=-50000 & genera<=1000000 &  constructed_net_worth>=-50000 & constructed_net_worth<=1000000 
drop genera


*CPI adjustment
foreach var in wealth constructed_net_worth family_li {
	gen real_`var'=.
	replace real_`var'=`var'*1 if year==1997
	replace real_`var'=`var'/166.600*160.500 if year==1999
	replace real_`var'=`var'/177.100*160.500 if year==2001
	replace real_`var'=`var'/183.960*160.500 if year==2003
	replace real_`var'=`var'/195.300*160.500 if year==2005
	replace real_`var'=`var'/207.342*160.500 if year==2007
	replace real_`var'=`var'/214.537*160.500 if year==2009
	replace real_`var'=`var'/224.939*160.500 if year==2011
	replace real_`var'=`var'/232.957*160.500 if year==2013
	replace real_`var'=`var'/237.017*160.500 if year==2015
	replace real_`var'=`var'/245.120*160.500 if year==2017
	replace real_`var'=`var'/255.7*160.500 if year==2019
	
	
	drop `var'

	}

rename real_wealth wealth
rename real_constructed_net_worth constructed_net_worth
rename real_family_li family_li

*As in EJ
replace constructed_net_worth=constructed_net_worth/10000
replace wealth=wealth/10000
replace family_li=family_li/10000

*Use the PSID variable after 2009
replace wealth=constructed_net_worth if year>=2009

*Drop tails...
sum wealth if hea==10,detail
keep if wealth<r(p99) & wealth>r(p1)
********************************************************************************
*Create relationship variables
********************************************************************************
*Create variable for spouse and unmarried partner
gen spouse=0
replace spouse=1 if hea==20
bysort id1968 year:egen mspouse=max(spouse)


gen partner=0
replace partner=1 if hea==22  | hea==88
bysort id1968 year:egen mpartner=max(partner)


*Keep only household heads
*keep if hea==10

*Generate lagged variables 
bysort id1968 pernum (year): gen pwealth=wealth[_n-1] if year==year[_n-1]+2
bysort id1968 pernum (year): gen pwealth2=wealth[_n-2] if year==year[_n-2]+4
bysort id1968 pernum (year): gen pms=ms_head[_n-1] if year==year[_n-1]+2
bysort id1968 pernum (year): gen ppartner=mpartner[_n-1] if year==year[_n-1]+2
bysort id1968 pernum (year): gen pspouse=mspouse[_n-1] if year==year[_n-1]+2
bysort id1968 pernum (year): gen phea=hea[_n-1] if year==year[_n-1]+2


label variable family_li "Current total household income"
label variable pwealth "Couple's net worth before divorce/breakup"

*Regressions
eststo clear
eststo:qreg wealth family_li pwealth if mr_status_change==3  & ms_head>=4 & ms_head<=4 & pspouse==1 & hea==10 
       su wealth if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	
eststo:qreg wealth family_li pwealth if mr_status_change==3  & ms_head>=5 & ms_head<=5 & ppartner==1 & hea==10 
       su wealth if e(sample), mean 
       loc mymean: di %8.3f r(mean) 
       estadd loc mD `mymean', replace	

	   *Table with results
  #delimit ;
  esttab using "C:\Users\Fabio\Dropbox\py_mar_3\Analysis\output\reg_wealth.tex", label replace 
    b(4) se(4) star(* 0.10 ** 0.05 *** 0.01) stats(mD N, labels("Dependent variable mean" "Observations") fmt(%9.0f) ) eqlabels(None) nonotes depvars nomtitles nonum
  
  prehead(\begin{table}[H]\centering 
				\scriptsize
				\caption{Changes in net worth upon divorce and breakup}  
				\label{tab:assdiv}
				\resizebox{0.8\textwidth}{!}{
				\begin{tabular}{l*{@M}{c}} \toprule
			    \textit{Dependent variable:}&\multicolumn{2}{c}{Inividual net worth after divorce/breakup}\\
				\textit{Sample:}& Divorces & Breakups  \\)	
			 posthead(\midrule)
			postfoot(
				\bottomrule \noalign{\smallskip}
				\end{tabular}
				}
				\begin{minipage}{\textwidth}
				\scriptsize\smallskip "Results from median regressions, a particular instance of the median regression, using data from the 1999-2019 PSID. The sample includes household head (men or women) who experienced a divorce or a breakup in the two years before the interview. We regress their individual net worth after the breakup on the couple's net worth before the breakup. Net worth and income are expressed in 1997 10,000 dollars. Coefficients that are significantly different from zero are denoted by *10\%, **5\%, and ***1\%."  \\ 
				\end{minipage}
				\end{table}
				)
				  ;
	#delimit cr
	

binscatter wealth pwealth pwealth if pwealth<=100 & pwealth>=0 & mr_status_change==3 & hea==10 & pspouse==1 & ms_head==5,median line(none)
binscatter wealth pwealth pwealth if pwealth<=100 & pwealth>=0 & mr_status_change==1 & hea==10 & pspouse==1 & ms_head==1,median line(none)
binscatter wealth pwealth pwealth if pwealth<=100 & pwealth>=0 & mr_status_change==3 & hea==10 & (phea==10 | pspouse==1 | ppartner==1), line(none)

gen married=0
replace married=1 if ms_head==1

***************************************
*Event Studies - divorce
***************************************

gen eyear=0
replace eyear=year if  mr_status_change==3  & ms_head>=4 & ms_head<=4  & pspouse==1 & hea==10
bysort id1968 pernum (year):egen meyear=max(eyear)

gen time_event=year-meyear

*welath
binscatter wealth time_event if time_event>=-8 & time_event<=6 , line(noline) by(sex) median

*marital status
binscatter married time_event if time_event>=-8 & time_event<=6, line(noline) by(sex)


did_imputation wealth id  year meyear if time_event<10000   ,fe(id year age)  autosample  allhorizons pre(8) nose
#delimit ;
event_plot , default_look graph_opt(xtitle("Years since the event") ytitle("Coefficients") legend(off) 
 xlabel(-6 "(-6,-5)" -4 "(-4,-3)" -2 "(-2,-1)"
         0 "(0,1)"      2 "(2,3)"     4 "(4,5)"    6 "(6,7)" , 
		 format(%13.0fc) labsize(4pt)))  trimlead(7) trimlag(8);
#delimit cr

gen id= id1968*10000000000000+pernum

gen time_event100=time_event+100

reg wealth i.age i.state i.year i.edu  ib98.time_event100 if time_event>=-6 & time_event<=6

drop eyear meyear time_event time_event100
***************************************
*Event Studies - breakup
***************************************
gen eyear=0
replace eyear=year if  mr_status_change==3  & ms_head>=5 & ms_head<=5  & ppartner==1 & hea==10
bysort id1968 pernum (year):egen meyear=max(eyear)

gen time_event=year-meyear

*welath
binscatter wealth time_event if time_event>=-8 & time_event<=6, line(noline) by(sex) median

*marital status
binscatter married time_event if time_event>=-8 & time_event<=6, line(noline) by(sex)

gen time_event100=time_event+100

reg wealth i.age i.state i.year i.edu  ib98.time_event100 if time_event>=-6 & time_event<=6

drop eyear meyear time_event time_event100

