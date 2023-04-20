*************************
*Create history of cohabitation
*******************************

*Get individual+family data
clear all
import delimited "C:\Users\Fabio\Dropbox\JMP\empirical analysis\Analysis\assets.csv", asdouble 

*drop if relationhead==0

rename pernum ER30002
rename id1968 ER30001

*Modify relationship to head to merge 1st year and 1yr+ cohabitor
replace relationhead=22 if relationhead==88

*account for swaps
gen head=-50000
replace head=ER30002 if relationhead==10 & sequence_number==1
bysort ER30001 interview year:egen headp=max(head)
drop head
drop if headp<0

gen partn=-50000
replace partn=ER30002 if relationhead==22
bysort ER30001 interview year:egen partnp=max(partn)
drop partn

*now couple is created
gen key1=string(max(headp,partnp))+string(min(headp,partnp)) if headp>0 & partnp>0
destring(key1),gen(key)


gen key1a=string(max(headp,partnp)) if headp>0 & partnp>0
gen key1b=string(min(headp,partnp)) if headp>0 & partnp>0
destring(key1a),gen(keya)
destring(key1b),gen(keyb)

drop key1 key1a key1b


*Get if marriage
gen married=0
replace married=1 if relationhead==20

bysort ER30001 key:egen marriedp=max(married) if key!=.

*Drop cohabitation that ended in marriage
*keep if relationhead==10 & sequence_number==1
drop if marriedp==1
drop married marriedp

*Now get cohabitation end etc.
gen cohabit=0
replace cohabit=1 if key!=.
bysort ER30001 key:egen cohabitt=max(cohabit)

gen separate=0
sort ER30001 ER30002 year
*replace separate=1 if (keya[_n-1]==ER30002[_n-1] | keyb[_n-1]==ER30002[_n-1])& (key[_n-1]!=key[_n]) & ER30002[_n-1]==ER30002[_n] & key[_n-1]!=.
replace separate=1 if (keya[_n]==ER30002[_n] | keyb[_n]==ER30002[_n])& (key[_n]!=key[_n+1]) & ER30002[_n]==ER30002[_n+1] & key[_n]!=.
bysort ER30001 key:egen sepp=max(separate)
drop if cohabitt==0
drop if sepp==0

bysort ER30001 key:egen beg=min(year)
bysort ER30001 key:egen sep=max(year)
replace sep=sep+2

bysort ER30001 key:keep if _n==1

keep ER30001 keya keyb beg sep

gen ER30001sp=ER30001
rename keya ER30002
rename keyb ER30002sp
rename beg mary
rename sep divy

merge m:1 ER30001 ER30002 using D:\blasutto\Data\PSID\family_history\temp.dta

*Keep only if matched
keep if _merge==3

*Rehsape
keep wealth* assets* divy mary    wgt* age* wmovein* sex* relationhead* iwy* iwm* n_fu* ch_fu*  non_fu* ER30001 ER30002 ER30001sp ER30002sp
*gen id=_n
reshape long wealth assets age sex wmovein iwy iwm n_fu ch_fu wgt relationhead non_fu,i(ER30001 ER30002 ER30001sp ER30002sp) j(year)

**************
*Gen event

drop if age>65
gen event=year-divy
replace event=. if iwy<mary 


*Before saving adjust the cpi
**** Convert ivaraibles into real values ******
foreach var in wealth assets {
	gen real_`var'=.
	replace real_`var'=`var'*1 if year==1997
	replace real_`var'=`var'*166.600/160.500 if year==1999
	replace real_`var'=`var'*177.100/160.500 if year==2001
	replace real_`var'=`var'*183.960/160.500 if year==2003
	replace real_`var'=`var'*195.300/160.500 if year==2005
	replace real_`var'=`var'*207.342/160.500 if year==2007
	replace real_`var'=`var'*214.537/160.500 if year==2009
	replace real_`var'=`var'*224.939/160.500 if year==2011
	replace real_`var'=`var'*232.957/160.500 if year==2013
	replace real_`var'=`var'*237.017/160.500 if year==2015
	replace real_`var'=`var'*245.120/160.500 if year==2017
	
	
	drop `var'
	*rename real_`var' head_li_
	}


rename real_assets assets
rename real_wealth wealth

*some regularity
bysort ER30001 ER30002 ER30001sp ER30002sp:egen min=min(event)
bysort ER30001 ER30002 ER30001sp ER30002sp:egen max=max(event)

*Conditions for keep
bysort ER30001 ER30002 ER30001sp ER30002sp:egen minw=min(wealth)
bysort ER30001 ER30002 ER30001sp ER30002sp:egen maxw=max(wealth)

keep if maxw<811006.6 

*Keep if together
keep if year-mary>=0

*Keep if transition
keep if min<0
keep if max>=0


*Get adults in the household
gen netto=n_fu-ch_fu
gen net=1
replace net=0 if n_fu-ch_fu>=2 & event>=0
replace net=0 if n_fu-ch_fu>=3 & event<0
replace net=0 if netto==1 & event<0

keep if net==1


replace event=7 if event>=7
replace event=-7 if event<=-7
gen event1=event+1000

gen lw=log(wealth+100000)
reg wealth i.event1 i.year i.mary i.age if event>=-3 & event<=6
predict xb,xb
gen lwh=exp(xb)-100000


*Keep only if transition observed

egen id=group (ER30001 ER30002 ER30001sp ER30002sp)

keep if event>=-7 & event<=7


***Regression below
reg wealth ib999.event1 i.year i.mary i.age 


gen wealthl=-100000000
gen wealthm=-100000000


gen below=0
replace below=1 if event<0


bysort id:egen wealthpr1=mean(wealth) if below==1


egen mwe=pctile(wealthpr1) if below==1,p(75) 
replace mwe=-1000000 if mwe==.7
bysort id:egen mwe1=max(mwe)

replace wealthpr1=-1000000 if wealthpr1==.
bysort id:egen wealthpr1f=max(wealthpr1)
gen less1=0
replace less1=1 if wealthpr1f>mwe1



*keep and save to have it in python (and to nice graphs)

keep if event1>=994 & event1<=1006
gen wealth2=wealth
replace wealth2=wealth/netto if event>=0
replace wealth2=wealth/netto*2 if event<0
gen dif1=year-mary

keep wealth wealth2 assets event less1 year mary age wgt sex id dif1 net

gen event1=event+1000
gen event2=event1
replace event2=1000 if event1>=1000

reg wealth2 i.mary i.age i.year b999.event2    if less==1 

saveold sepe.dta,replace
