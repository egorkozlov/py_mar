*****************************************************
*Use PSID data to build estimate men income process
****************************************************

clear
set more off

import delimited C:\Users\Fabio\Dropbox\py_mar_3\Analysis\income_estimation\income_men_ch1_old.csv


********************************
*Destring+Clean some earnings variables
********************************
quietly destring *, replace ignore("NA")

gen miy=1968
gen may=1993

keep if year<=may & year>=miy

bysort pid (year): egen miny=min(year)
bysort pid (year): egen mage=min(age)
gen cohort=year-age


*********************************
*When us a variable topcoded?
*Below you find a list:
*******************************

*Head
gen topc=0
replace topc=65491 if year==1968 
replace topc=80001 if year==1969 
replace topc=99999 if year>=1970 & year<=1982 
replace topc=999999 if year>=1983 & year<=1992
replace topc=9999999 if year==1993 
replace topc=9999999 if year==1996
replace topc=999999 if year==1997 
replace topc=9999999 if year>=1999 & year<=2009 
replace topc=9999998 if year>=2011 



*Fit pareto distribution of non topcoded data, then predict the topcoded values
gen old_head_li=head_li
foreach i of numlist 1968/1993{
	
	*Get to 1%
	 quietly summarize old_head_li if year==`i' & old_head_li>0 & old_head_li<topc,detail
	 quietly g top=r(p99)
	*Fit pareto distribution
	 quietly paretofit old_head_li if year==`i' & old_head_li>0 & old_head_li<topc,x0(top)
	
	 *Change values if topcoded
	 quietly replace head_li=e(ba)/(e(ba)-1)*topc if old_head_li>=topc & year==`i' 
	 quietly drop top
}

********************************
*Clean some more variables
********************************

*Indicator whether is single
gen married=0
replace married=1 if ms==1

*Resize age to be consistent with the model
*and that we have info about previous year
replace age=age-21

*Rename person
rename pid person

*** Code the minumum wage *****
gen mwage=.
replace mwage=1.60 if year>=1968 & year<1974
replace mwage=2.00 if year>=1974 & year<1975
replace mwage=2.10 if year>=1975 & year<1976
replace mwage=2.30 if year>=1976 & year<1978
replace mwage=2.65 if year>=1978 & year<1979
replace mwage=2.90 if year>=1979 & year<1980
replace mwage=3.10 if year>=1980 & year<1981
replace mwage=3.35 if year>=1981 & year<1990
replace mwage=3.80 if year>=1990 & year<1991
replace mwage=4.25 if year>=1991 & year<1996
replace mwage=4.75 if year==1996
replace mwage=5.15 if year>=1997 & year<2006
replace mwage=5.85 if year==2007
replace mwage=6.55 if year==2008
replace mwage=7.25 if year>2008

**** Convert income into real values ******
foreach var in head_li mwage {
	gen real_`var'=.
	replace real_`var'=`var'*3.739003 if year==1968
	replace real_`var'=`var'*3.571429 if year==1969
	replace real_`var'=`var'*3.364116 if year==1970
	replace real_`var'=`var'*3.195488 if year==1971
	replace real_`var'=`var'*3.09466 if year==1972
	replace real_`var'=`var'*2.985949 if year==1973
	replace real_`var'=`var'*2.724359 if year==1974
	replace real_`var'=`var'*2.437859 if year==1975
	replace real_`var'=`var'*2.284946 if year==1976
	replace real_`var'=`var'*2.172061 if year==1977
	replace real_`var'=`var'*2.033493 if year==1978
	replace real_`var'=`var'*1.861314 if year==1979
	replace real_`var'=`var'*1.634615 if year==1980
	replace real_`var'=`var'*1.462156 if year==1981
	replace real_`var'=`var'*1.350636 if year==1982
	replace real_`var'=`var'*1.302349 if year==1983
	replace real_`var'=`var'*1.248776 if year==1984
	replace real_`var'=`var'*1.206244 if year==1985
	replace real_`var'=`var'* 1.2560146 if year==1986
	replace real_`var'=`var'* 1.2544524 if year==1987
	replace real_`var'=`var'*1.099138 if year==1988
	replace real_`var'=`var'*1.05198 if year==1989
	replace real_`var'=`var'*1 if year==1990
	replace real_`var'=`var'*0.9465479 if year==1991
	replace real_`var'=`var'*0.9219089 if year==1992
	replace real_`var'=`var'*0.8928571 if year==1993
	
	*Adjust to make comparable to 1997
	replace real_`var'=real_`var'*1.228
	
	
	drop `var'
	}


g li=real_head_li
gen ln_ly=ln(real_head_li)


*Normalization
replace ln_ly=log(exp(ln_ly)/2000)-2.362130
*************************************************************
*Follow Heathcote, Perri and Violante (2010) to subset data
*************************************************************


*Following Daly et al. 2022: 
*we drop the first three and last three earnings (or wage) observations if an 
*individual's first record is after 1979 (1968 for us) and the last record is prior to 1992, 
*as well as the three earnings (or wage) observations before and after missing earnings (or wage) records.

*First and last observations...
bysort person:egen minny=min(year)
bysort person:egen maxxy=max(year)


gen agei=age+100
reg ln_ly i.agei i.year i.state
predict residuals,residuals

gen diffe_b=1000000
gen diffe_e=1000000
replace diffe_b=year-minny if minny>miy
replace diffe_e=year-maxxy+100 if maxxy<may

*...missing observations
drop if ln_ly==.
gen missing=0
bysort person:replace missing=1 if (year[_n]-year[_n-1]>1) | (year[_n+1]-year[_n]>1)  & year>minny & year<maxxy

bysort person (year):gen min=sum(missing)
bysort person min:egen minny_m=min(year)
bysort person min:egen maxxy_m=max(year)
gen diffe_b_m=1000000
gen diffe_e_m=1000000
replace diffe_b_m=year-minny_m if minny!=minny_m
replace diffe_e_m=year-maxxy_m+100 if maxxy!=maxxy_m

*Drop observations here
drop if diffe_b>=0 & diffe_b<=2
drop if diffe_e>=100 & diffe_e<=100
drop if diffe_b_m>=0 & diffe_b_m<=2
drop if diffe_e_m>=98 & diffe_e_m<=100


*Consider just men that are working
drop if real_head_li==0 | head_hrs==0

*Keep only between age 20 to age 64
keep if age>=0 & age<=44

*No less tha half the minimum wage
drop if li<0.5*real_mwage*head_hrs

*Drop if less than 260 hours a year
drop if head_hrs<260 & real_head_li>0

*Drop the bottom 0.5% of earnings
centile(ln_ly), centile(0.5)
local cent=r(c_1)
drop if ( ln_ly < `cent' )

replace edu=-1 if edu==9 | edu==.
bysort person (year):egen medu=max(edu)
*drop if medu>=7 & age<=5

*Very good...with 150 hours or 200, 300 hourss
*Down: the more you trim, the more precise it gets
bysort person (year):gen dd_ln_ly=ln_ly-ln_ly[_n-1] if year==year[_n-1]+1
bysort person (year):gen di_ln_ly=real_head_li/real_head_li[_n-1] if year==year[_n-1]+1
drop if (di_ln_ly> 3.5   | di_ln_ly< .3333333333) & di_ln_ly!=.


**************************************************************************
****** ESTIMATE VARIANCE OF WAGE PROCESS WITH SELECTION INTO WORK  *******
**************************************************************************

*Variables for regression
gen age2=age^2
gen age3=age^3

tab(state), gen(std)
tab year, gen(yrd)

gen newc=0
replace newc=1 if youngest_child>=0 & youngest_child<=3
replace ch_fu=3 if ch_fu>=3

reg ln_ly age age2   i.year i.state i.ms newc i.ch_fu 


*Store trends
gen t1=_b[age]
gen t2=_b[age2]
gen bsingle=_b[2.ms]
*gen t3=_b[age3]


bys person (year): gen  d_ch_fu= ch_fu-ch_fu[_n-1]+10 if year==year[_n-1]+1
bys person (year): gen  d_newc= newc-newc[_n-1]+10 if year==year[_n-1]+1


gen edu2=edu+100

bys person (year): gen  d_edu= edu2-edu2[_n-1]+10 if year==year[_n-1]+1
bys person (year): gen  d_ms= ms-ms[_n-1]+10 if year==year[_n-1]+1

reg dd_ln_ly i.age i.year i.ms i.state newc i.ch_fu i.d_ch_fu i.edu2  



*Generate residuals for GMM
predict y_m if e(sample)==1, resid
predict y_m_hat if e(sample)==1, xb

sort person year
gen y_m2= (y_m)^2
sort person year
gen u_0=y_m 
bys person: gen  u_1= y_m[_n-1] if year==year[_n-1]+1 
bys person: gen  u_2= y_m[_n-2] if year==year[_n-2]+2 
bys person: gen  u_3= y_m[_n-3] if year==year[_n-3]+3 
bys person: gen  u_4= y_m[_n-4] if year==year[_n-4]+4
bys person: gen  u_5= y_m[_n-5] if year==year[_n-5]+5 
************************************
*Get Variance shocks with residuals
**************************************
gen expect=u_2*(u_0+u_1+u_2+u_3+u_4) 

*Get SE and save result
*tabstat expect  , stat(semean)
sum expect
g var_per_t=r(mean)



*keep if cohort>=1980 & cohort<=1984
*Implied inequality over age
bysort age:egen vagevae=sd(ln_ly) 
forvalues i=0(1)44{
	
	quietly sum ln_ly if age==`i' ,detail
	replace vagevae=r(Var) if age==`i'
}

nlsur (vagevae={v_int}+age*var_per_t) if var_per_t!=-. & age<=15 & age>=5
gen v_int=_b[/v_int]
g var_age=v_int+(age)*var_per

*Create needed variables
egen wmean=mean(ln_ly),by(age)
egen wvar=sd(ln_ly),by(age)
replace wvar=wvar^2

gen agem=age+20
#delimit ;
twoway (line var_age agem,lwidth(1) color(black)) (scatter wvar agem, color(red) msymbol(Oh) msize(2.5))  if age>=5, 
xtitle("Age") ytitle("variance") xlabel(25(10)65)  ylabel(0.0(0.2)1.2) 
		legend(order(1 "Estimated" 2 "Data") region(style(none)))
		graphregion(color(white)) bgcolor(white)
;
#delimit cr
graph export "C:\Users\Fabio\Dropbox\py_mar_3\Analysis\output\variance_age.png",replace

*Check how earnings profile look like
bysort age: egen earn_age=mean(ln_ly)
nlsur (ln_ly={vint}+(age)*t1+(age2)*t2) if age<=15 & age>=5
gen t0m=_b[/vint]
gen pred=t0m+(age)*t1+(age2)*t2
twoway (scatter earn_age age) (scatter pred age)


*Save stuff
file open myfile using "C:\Users\Fabio\Dropbox\py_mar_3\Analysis\output\earnings_params_m.csv", write replace

foreach v of varlist t0m t1 t2 v_int expect bsingle {
       sum `v'
	   local vmean: display %10.8fc `r(mean)'
	   file write myfile "`vmean' " _n
}


*Also with housework
sum total_hsw if ms==2 ,detail
local vmean: display %10.8fc `r(mean)'
file write myfile "`vmean' " _n
file close myfile
*****************************************************
*Use PSID data to build estimate male women process
****************************************************

*****************************************
*A Note: if we consider houry wage (tot labor incom/hours year) then
*the persistene shock for men and women is similar, otherwise if you
*consider the whole year labor income, the the shock results higher for
*women

clear all
*set more off

import delimited C:\Users\Fabio\Dropbox\py_mar_3\Analysis\income_estimation\income_women_ch1_old.csv

********************************
*Destring+Clean some earnings variables
********************************
quietly destring *, replace ignore("NA")

bysort pid (year): egen miny=min(year)
bysort pid (year): egen mage=min(age)
gen cohort=year-age

*College or noncollege?
gen miy=1968
gen may=1993

keep if year<=may & year>=miy
*********************************
*When us a variable topcoded?
*Below you find a list:
*******************************

*Head
gen topc=0
replace topc=65491 if year==1968 
replace topc=80001 if year==1969 
replace topc=99999 if year>=1970 & year<=1982 
replace topc=999999 if year>=1983 & year<=1992
replace topc=9999999 if year==1993 
replace topc=9999999 if year==1996
replace topc=999999 if year==1997 
replace topc=9999998 if year>=2011 

*Partner 
gen topcf=0
replace topcf=15001 if year==1968 
replace topcf=12459 if year==1969 
replace topcf=99999 if year>=1970 & year<=1983 
replace topcf=999999 if year>=1984 & year<=1992 
replace topcf=9999999 if year==1993 
replace topcf=9999999 if year==1996 
replace topcf=999999 if year==1997 
replace topcf=9999998 if year==1999 
replace topcf=9999999 if year>=2001 & year<=2009 
replace topcf=99999998 if year>=2011 



	 *Fit pareto distribution of non topcoded data, then predict the topcoded values
     gen old_head_li=head_li
	 gen old_partner_li=partner_li
     foreach i of numlist 1968/1993{
	
	*Get to 1%
	 quietly summarize old_head_li if year==`i' & old_head_li>0 & old_head_li<topc,detail
	 quietly g top=r(p99)
	*Fit pareto distribution
	 quietly paretofit old_head_li if year==`i' & old_head_li>0 & old_head_li<topc,x0(top)
	
	 *Change values if topcoded
	 quietly replace head_li=e(ba)/(e(ba)-1)*topc if old_head_li>=topc & year==`i' 
	 quietly drop top
	 
	 
	 	*Get to 1%
	 quietly summarize old_partner_li if year==`i' & old_partner_li>0 & old_partner_li<topcf,detail
	 quietly g top=r(p99)
	*Fit pareto distribution
	 quietly paretofit old_partner_li if year==`i' & old_partner_li>0 & old_partner_li<topcf,x0(top)
	
	 *Change values if topcoded
	 quietly replace partner_li=e(ba)/(e(ba)-1)*topcf if old_partner_li>=topcf & year==`i' 
	 quietly drop top
}



********************************************************
*Take the right observation for age, income and hours
********************************************************
gen wage=.
replace wage=head_li if hea==1 & year<1983 & sex==2
replace wage=head_li if hea==10 & year>=1983 & sex==2
replace wage=partner_li if hea==1 & year<1983 & sex==1
replace wage=partner_li if hea==10 & year>=1983 & sex==1

gen age_h=age

gen agew=.
replace agew=age if hea==1 & year<1983 & sex==2
replace agew=age if hea==10 & year>=1983 & sex==2
replace agew=age_w if hea==1 & year<1983 & sex==1
replace agew=age_w if hea==10 & year>=1983 & sex==1

gen hours=.
replace hours=head_hrs if hea==1 & year<1983 & sex==2
replace hours=head_hrs if hea==10 & year>=1983 & sex==2
replace hours=part_hrs if hea==1 & year<1983 & sex==1
replace hours=part_hrs if hea==10 & year>=1983 & sex==1

gen educ=.
replace educ=edu if hea==1 & year<1983 & sex==2
replace educ=edu if hea==10 & year>=1983 & sex==2
replace educ=eduf if hea==1 & year<1983 & sex==1 & eduf>0
replace educ=eduf if hea==10 & year>=1983 & sex==1 & eduf>0

gen husband_li=head_li if hea==1 & year<1983 & sex==1
replace husband_li=head_li if hea==10 & year>=1983 & sex==1
gen husband_edu=edu if hea==1 & year<1983 & sex==1
replace husband_edu=edu if hea==10 & year>=1983 & sex==1
gen husband_college=.
replace husband_college=0 if husband_edu<16 
replace husband_college=1 if husband_edu>=16 & husband_edu<=30

**********************************
*Some more clearning
*********************************


*Rename person
rename pid person



*Indicator whether is single
gen married=0
replace married=1 if ms==1
drop if married==0 & sex==1

*Rename age
replace age=agew
* if married==1


*Resize age to be consistent with the model
*and that we have info about previous year
replace age=age-21



*** Code the minumum wage *****
gen mwage=.
replace mwage=1.60 if year>=1968 & year<1974
replace mwage=2.00 if year>=1974 & year<1975
replace mwage=2.10 if year>=1975 & year<1976
replace mwage=2.30 if year>=1976 & year<1978
replace mwage=2.65 if year>=1978 & year<1979
replace mwage=2.90 if year>=1979 & year<1980
replace mwage=3.10 if year>=1980 & year<1981
replace mwage=3.35 if year>=1981 & year<1990
replace mwage=3.80 if year>=1990 & year<1991
replace mwage=4.25 if year>=1991 & year<1996
replace mwage=4.75 if year==1996
replace mwage=5.15 if year>=1997 & year<2006
replace mwage=5.85 if year==2007
replace mwage=6.55 if year==2008
replace mwage=7.25 if year>2008

**** Convert income into real values ******
foreach var in  wage husband_li mwage {
	gen real_`var'=.
	replace real_`var'=`var'*3.739003 if year==1968
	replace real_`var'=`var'*3.571429 if year==1969
	replace real_`var'=`var'*3.364116 if year==1970
	replace real_`var'=`var'*3.195488 if year==1971
	replace real_`var'=`var'*3.09466 if year==1972
	replace real_`var'=`var'*2.985949 if year==1973
	replace real_`var'=`var'*2.724359 if year==1974
	replace real_`var'=`var'*2.437859 if year==1975
	replace real_`var'=`var'*2.284946 if year==1976
	replace real_`var'=`var'*2.172061 if year==1977
	replace real_`var'=`var'*2.033493 if year==1978
	replace real_`var'=`var'*1.861314 if year==1979
	replace real_`var'=`var'*1.634615 if year==1980
	replace real_`var'=`var'*1.462156 if year==1981
	replace real_`var'=`var'*1.350636 if year==1982
	replace real_`var'=`var'*1.302349 if year==1983
	replace real_`var'=`var'*1.248776 if year==1984
	replace real_`var'=`var'*1.206244 if year==1985
	replace real_`var'=`var'* 1.2560146 if year==1986
	replace real_`var'=`var'* 1.2544524 if year==1987
	replace real_`var'=`var'*1.099138 if year==1988
	replace real_`var'=`var'*1.05198 if year==1989
	replace real_`var'=`var'*1 if year==1990
	replace real_`var'=`var'*0.9465479 if year==1991
	replace real_`var'=`var'*0.9219089 if year==1992
	replace real_`var'=`var'*0.8928571 if year==1993
	
	*Adjust to make comparable to 1997
	replace real_`var'=real_`var'*1.228
	
	
	drop `var'
	*rename real_`var' head_li_
	}

gen ln_ly=ln(real_wage)
g li=real_wage

*Normalization
replace ln_ly=log(exp(ln_ly)/2000)-2.362130

*************************************************************
*Follow Heathcote, Perri and Violante (2010) to subset data
*************************************************************



*Following Daly et al. 2022: 
*we drop the first three and last three earnings (or wage) observations if an 
*individual's first record is after 1979 (1968 for us) and the last record is prior to 1992, 
*as well as the three earnings (or wage) observations before and after missing earnings (or wage) records.

*First and last observations...
bysort person:egen minny=min(year)
bysort person:egen maxxy=max(year)


gen agei=age+100
reg ln_ly i.agei i.year i.state
predict residuals,residuals

gen diffe_b=1000000
gen diffe_e=1000000
replace diffe_b=year-minny if minny>miy
replace diffe_e=year-maxxy+100 if maxxy<may

*...missing observations
gen missing=0
bysort person:replace missing=1 if (year[_n]-year[_n-1]>1) | (year[_n+1]-year[_n]>1)  & year>minny & year<maxxy

bysort person (year):gen min=sum(missing)
bysort person min:egen minny_m=min(year)
bysort person min:egen maxxy_m=max(year)
gen diffe_b_m=1000000
gen diffe_e_m=1000000
replace diffe_b_m=year-minny_m if minny!=minny_m
replace diffe_e_m=year-maxxy_m+100 if maxxy!=maxxy_m

*Drop observations here
drop if diffe_b>=0 & diffe_b<=2
drop if diffe_e>=100 & diffe_e<=100
drop if diffe_b_m>=0 & diffe_b_m<=2
drop if diffe_e_m>=98 & diffe_e_m<=100


*Keep only certain ages
keep if age>=-2 & age<=44

*Drop if inconsistency (ex zero hours and positive wage and vice versa)
drop if hours==0 & real_wage>0
drop if hours>0 & real_wage==0

*No less tha half the minimum wage
drop if li<0.5*real_mwage*hours & hours>0


*Drop the bottom 0.5% of earnings
centile(ln_ly), centile(0.5)
drop if ( ln_ly < r(c_1) )

*Drop if wage changes are too crazy 
bysort person (year):gen d_ln_ly=ln_ly-ln_ly[_n-1] if year==year[_n-1]+1
bysort person (year):gen di_ln_ly=real_wage/real_wage[_n-1] if year==year[_n-1]+1 
drop if (di_ln_ly> 3.5  | di_ln_ly< .3333333) & di_ln_ly!=. 
 
*Limit on hrs 
*drop if hours<=100 & real_wage>0 
drop if hours<=260 & ln_ly!=.


replace edu=-1 if edu==9 | edu==.
bysort person (year):egen medu=max(edu)
*drop if medu>=7 & age<=5

*****Generate female labor force participation
gen fp=0
replace fp=1 if hours>0


*Age
gen age2=age^2

*Generate the isntrument
g mortgage1=0
replace mortgage1=1 if mort>=1 & mort<=2

g mortgage2=0
replace mortgage2=1 if mort==2

gen newc=0
replace newc=1 if youngest_child>=1 & youngest_child<=3
replace ch_fu=3 if ch_fu>=3

g agem=age+2

*get labor experience
gen missed=.
replace missed=0 if fp==1
replace missed=1 if fp==0
bysort person (year):gen sfp=sum(missed)+fp

sort person year

bysort person:gen obs=_n



*Probit for participation (need for different model for nocollege, otherwise ggap too large compared to data.)


probit fp  age age2 i.year i.state i.ms  newc i.ch_fu  mortgage1 

*labels for the graph
label var age2 "\$\iota^f_2\$"
label var age "\$\iota^f_1\$"
label var mortgage1 "Mortgage"



   #delimit ;
  esttab using "C:\Users\Fabio\Dropbox\py_mar_3\Analysis\output\prb_wom", label
  b(2) not  booktabs  nomtitles  replace style(tex) keep(mortgage1 age age2) nolines nonum 
  eqlabels(none)
  prefoot("Survey Year Fixed Effects  & \checkmark   \\
            State Fixed Effects  & \checkmark  \\
		  Demographic Controls  & \checkmark  \\
			\hline")
			postfoot(\hline)
  prehead("{\def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi}
             \begin{tabular}{l*{1}{c}}
			 \toprule
            &\multicolumn{1}{c}{(1)}       \\
            \midrule
            \textsc{Dep. Variable:} &  \\\textsc{Female Labor Force Participation} & \\ & \\");
	#delimit cr
	
	eststo clear


predict gm if e(sample)==1
replace gm=-gm
gen     lambdal=normalden(gm)/(1-normal(gm))

*Now the full regerssion
reg ln_ly age age2  i.year i.state i.ms newc i.ch_fu   lambdal
gen bsingle=_b[2.ms]
drop gm 
predict ln_ly_ha, xb 


*Store trend
gen t1=_b[age]
gen t2=_b[age2]


*Now the full regerssionv
bys person (year): gen  d_ch_fu= ch_fu-ch_fu[_n-1]+10 if year==year[_n-1]+1
bys person (year): gen  d_newc= newc-newc[_n-1]+10 if year==year[_n-1]+1
bys person (year): gen  d_obs= obs-obs[_n-1]+10 if year==year[_n-1]+1
bys person (year): gen  d_sfp= sfp-sfp[_n-1]+10 if year==year[_n-1]+1 
gen edu2=edu+100
gen age100=age+100

probit fp  i.age100 i.year i.ms i.state newc i.ch_fu i.d_ch_fu i.d_newc i.edu2 mortgage1 

*Keep only if in sample
drop if e(sample)==0 

 
*Get inverse of the mills ratio
predict gm if e(sample)==1
replace gm=-gm
gen     lambda=normalden(gm)/(1-normal(gm))
bysort person (year): gen d_lambda=lambda-lambda[_n-1] if year==year[_n-1]+1

reg d_ln_ly i.age100 i.year i.ms i.state newc i.ch_fu i.d_ch_fu i.edu2 d_lambda  



*Get residuals and predition
predict y_m if e(sample)==1, resid






gen u_0=y_m
bys person: gen  u_p2= y_m[_n+2] if year==year[_n+2]-2 
bys person: gen  u_p1= y_m[_n+1] if year==year[_n+1]-1  
bys person: gen  u_1= y_m[_n-1] if year==year[_n-1]+1  
bys person: gen  u_2= y_m[_n-2] if year==year[_n-2]+2 
bys person: gen  u_3= y_m[_n-3] if year==year[_n-3]+3  
bys person: gen  u_4= y_m[_n-4] if year==year[_n-4]+4  
bys person: gen  u_5= y_m[_n-5] if year==year[_n-5]+5 





*****************************************************************************
*Get the variables with GMM: this follows Heatcote, Perri Violante 2010, p25
****************************************************************************
**Create variables for earnings growth


********************************** 
*With Participation correction 
********************************* 
gen expect=u_2*(u_0+u_1+u_2+u_3+u_4)  
  
/*
*GMM 
#delimit ; 
 
gmm  (u_0 -{std_pers}*lambda) 
	 (expect- {var_pres}-{std_pers}^2*lambda*gm)
     , igmm quickd  
	winitial(identity); 
# delimit cr 
 
g var_pers=_b[/var_pres] 
*/
g var_pers=expect
*Implied inequality over age
bysort age:egen vagevae=sd(ln_ly)
forvalues i=0(1)44{
	
	quietly sum ln_ly if age==`i' ,detail
	replace vagevae=r(Var) if age==`i'
}
nlsur (vagevae={v_int}+age*var_pers)  if age<=15 & age>=5 & var_pers!=.
gen v_int=_b[/v_int]
g var_age=v_int+(age)*var_pers

*Create needed variables
egen wmean=mean(ln_ly),by(age)
egen wvar=sd(ln_ly),by(age)
replace wvar=wvar^2

twoway (scatter wvar age) (scatter var_age age)

*Check how earnings profile look like
bysort age: egen earn_age=mean(ln_ly)
nlsur (ln_ly_ha={vint}+(age)*t1+(age2)*t2) if age<=15 & age>=5
gen t0=_b[/vint]
gen pred=t0+(age)*t1+(age2)*t2
twoway (scatter earn_age age) (scatter pred age)


*Save stuff
file open myfile using "C:\Users\Fabio\Dropbox\py_mar_3\Analysis\output\earnings_params_f.csv", write replace

foreach v of varlist t0 t1 t2 v_int expect bsingle {
       sum `v'
	   local vmean: display %10.8fc `r(mean)'
	   file write myfile "`vmean' " _n
}

*Also with housework
sum total_hsw if ms==2,detail
local vmean: display %10.8fc `r(mean)'
file write myfile "`vmean' " _n




************************************ 
*Women and Husbands - works if college=2 
************************************ 
 
*Sample with husbands 
gen males=1 
replace males=0 if married==0 

gen ln_ly_hu=log(head_li)  
 
centile(ln_ly_hu), centile(0.5) 
replace males=0 if  ( ln_ly_hu < r(c_1) ) & married==1 
 
*Drop if wage changes are too crazy 
gen eln_ly=ln(li) 
bysort person (year):gen dii_ln_ly=head_li/head_li[_n-1] if year==year[_n-1]+1 
replace males=0 if (dii_ln_ly> 3.5  | dii_ln_ly< .333333) & di_ln_ly!=. & married==1 
 

 
replace ln_ly_hu=log(exp(ln_ly_hu)/2000)-2.1317884
bysort person (year): gen first=1 if _n==1 
 
*Correlation in earnings-all marries(too few people when we know they just got together) 
 
corr ln_ly ln_ly_hu if males==1 
 
*Share income women 
gen share=exp(ln_ly)/(exp(ln_ly_hu)+exp(ln_ly)) if males==1 
 
*Assortative mating by earnings 
gen colle=0 
replace colle=1 if educ>=16 
 
 
*Reg husband on wife
reg ln_ly_hu ln_ly   if males==1 &   first==1 & married==1 
predict mal,resid 
sum mal  if males==1 &   first==1 & married==1 
local vmean: display %10.8fc `r(sd)'
file write myfile "`vmean' " _n

gen bwife=e(b)[1,1]
sum bwife
local vmean: display %10.8fc `r(mean)'
file write myfile "`vmean' " _n

*Reg wife on husband
reg ln_ly ln_ly_hu    if males==1 & first==1 & married==1
predict fem,resid 
sum fem  if males==1 &  first==1 & married==1
local vmean: display %10.8fc `r(sd)'
file write myfile "`vmean' " _n

gen bhu=e(b)[1,1]
sum bhu
local vmean: display %10.8fc `r(mean)'
file write myfile "`vmean' " _n


file close myfile
 
 