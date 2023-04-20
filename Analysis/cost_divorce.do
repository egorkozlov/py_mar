***************************
*Cost of Divorce in the PSID
*****************************



*********************************
*Get assets
*********************************
clear all
import delimited "C:\Users\Fabio\Dropbox\JMP\empirical analysis\Analysis\assets.csv", asdouble 

global clean1 twh other_debt2 ba_a2 ba_d2 impw2_deb impw2_ass other_debt busval real_estate
global clean2 penval busval house real_estate other_debt stocks bonds cash carval mort1 mort2

foreach val in $clean1{ 

replace `val'="." if `val'=="NA"
destring `val',gen(`val'1)
drop `val'
rename `val'1 `val'

}


foreach val in $clean2{

replace `val'=. if `val'>=999999998
replace `val'=. if `val'<=-99999998
                          
}
replace mort1=. if mort1>=9999998
replace mort2=. if mort2>=9999998
replace stocks=. if stocks>=99999998

*Generate variations of wealth
g wealth1=twh
g wealth2=cash+bonds+stocks+busval+penval+house+real_estate+carval-other_debt-mort1-mort2
g wealth3=house+cash+bonds+stocks+penval+impw2_ass-impw2_deb+ba_a2-ba_d2+carval-other_debt2-mort1-mort2

gen wealth=wealth2
replace wealth=wealth3 if wealth==.

*gen wealth=wealth1
*replace wealth=wealth2 if wealth==.

g assets=cash+bonds+stocks+busval+penval+house+real_estate+carval
replace assets=house+cash+bonds+stocks+penval+impw2_ass-impw2_deb+ba_a2-ba_d2+carval if assets==.

keep pernum id1968 wealth assets year age sex wmovein iwy iwm n_fu ch_fu non_fu wgt relationhead 

*As in BPS
keep if wealth>=-1000000 & wealth<20000000


reshape wide wealth assets age sex wmovein iwy iwm n_fu ch_fu wgt relationhead non_fu, i(id1968 pernum) j(year)

rename pernum ER30002
rename id1968 ER30001

saveold D:\blasutto\Data\PSID\family_history\temp.dta,replace

clear all
cd "C:\Users\Fabio\Dropbox\JMP\empirical analysis\Analysis"


*********************************
*MANIPULATE MARITAL HISTORY FILE
*********************************

*first, import data on marriages ad keep what is needed
use D:\blasutto\Data\PSID\family_history\fhist.dta

*Keep if divorced
keep if MH12==4

*Keep if divorce happened 1999-2017
keep if MH14>=1996 & MH14<=2017

*Keep if now difference between divorce and separation
*keep if (MH14==MH16) | MH16>2017


*Generate variables of interest
gen divy=MH14
gen divm=MH13
gen mary=MH11
gen marm=MH10
gen sepm=MH15 if MH15<=12
gen sepy=MH16 if MH16<=2020
gen MH2sp=MH7
gen MH3sp=MH8

*Creat duplicates for spouse
expand 2
bysort MH2 MH3 MH9:gen sp=_n
replace MH2sp=MH2 if sp==2
replace MH3sp=MH3 if sp==2
replace MH2=MH7 if sp==2
replace MH3=MH8 if sp==2

*Duplicate again if main respondent accoring to number of divorces
*bysort MH2 MH3 sp:egen ndiv=count(MH2)
*replace ndiv=1 if sp==2

*Drop if spouse neve in the individual file
*drop if MH3>=800



keep MH2 MH3 MH2sp MH3sp sp divy divm mary marm sepy sepm  
rename MH2 ER30001
rename MH3 ER30002
rename MH2sp ER30001sp
rename MH3sp ER30002sp

*Eliminate some doubles
bysort ER30001 ER30002 ER30001sp ER30002sp:gen el=_n
drop if el>1
drop el

*********************************
*MERGE WITH INDIVIDUAL FILE
*********************************
merge m:1 ER30001 ER30002 using D:\blasutto\Data\PSID\family_history\temp.dta

*Keep only if matched
keep if _merge==3

*Rehsape
keep wealth* assets* divy divm mary sepy marm sepm sp wgt* age* wmovein* sex* relationhead* iwy* iwm* n_fu* ch_fu*  non_fu* ER30001 ER30002 ER30001sp ER30002sp
*gen id=_n
reshape long wealth assets age sex wmovein iwy iwm n_fu ch_fu wgt relationhead non_fu,i(ER30001 ER30002 ER30001sp ER30002sp) j(year)

**************
*Gen event
gen event_raw=(iwy-divy)*12+iwm-divm if divm<=12


drop if age>65
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

	}

rename real_assets assets
rename real_wealth wealth
gen event=year-divy
replace event=. if iwy*100+iwm<mary*100+marm & marm<=12
replace event=. if iwy<mary & marm>12

*some regularity
bysort ER30001 ER30002 ER30001sp ER30002sp:egen min=min(event)
bysort ER30001 ER30002 ER30001sp ER30002sp:egen max=max(event)

*Conditions for keep
bysort ER30001 ER30002 ER30001sp ER30002sp:egen minw=min(wealth)
bysort ER30001 ER30002 ER30001sp ER30002sp:egen maxw=max(wealth)


*Dont keep if too rich, 4 percent
keep if maxw< 817639.3


*3795160,1506276,1088860,817639.3,666594.8

*-74281.4,-19374.53
*888846.4,1174157,1607422
*4157072
*342340.4
*700308.8
*538000 
*4157072


*Get adults in the household
gen netto=n_fu-ch_fu
gen net=1
replace net=0 if n_fu-ch_fu>=2 & event>=0
replace net=0 if n_fu-ch_fu>=3 & event<0
replace net=0 if netto==1 & event<0


*Gen id1
egen id=group (ER30001 ER30002 ER30001sp ER30002sp)

*Keep if one person lost after divorce
keep if net==1




drop if event>6
drop if event<-6
replace event=7 if event>=7
replace event=-7 if event<=-7
gen event1=event+1000

*Keep if before and after
bysort id:egen min1=min(event)
bysort id:egen max1=max(event)
drop if min1>=0
drop if max1<0



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




keep wealth assets event less1 year mary age wgt sex id  net

gen event1=event+1000

gen wealth2=wealth
replace wealth2=wealth*2 if event>=0

gen event2=event1
replace event2=1000 if event1>=1000

reg wealth2 i.mary i.age i.year b999.event2    if less==1 

saveold dive.dta,replace

