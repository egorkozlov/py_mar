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
gen length=COHABINT

*generate outcome
gen sep=0
replace sep=1 if COHOUT==4

gen mar=0
replace mar=1 if COHOUT==2 | COHOUT==3

*keep only if at least one month spell
drop if COHABINT<=1

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
replace unil=0 if state=="Missouri" 
replace unil=1973 if state=="Montana" 
replace unil=1972 if state=="Nebraska" 
replace unil=1967 if state=="Nevada" 
replace unil=1971 if state=="New Hampshire" 
replace unil=0 if state=="New Jersey" 
replace unil=1933 if state=="New Mexico" 
replace unil=0 if state=="New York" 
replace unil=0 if state=="North Carolina" 
replace unil=1971 if state=="North Dakota" 
replace unil=0 if state=="Ohio" 
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

*education
gen coll=0
replace coll=1 if EDUCAT>=16

*recode to make it comparable with dates
replace unil=(unil-1900)*12+1 if unil!=0
replace eq=(eq-1900)*12+1 if eq!=0

****************************************************
*****************************************************
**Preparte to Expand the data***
ge id=_n
*replace length=length+1 if mar==1
*replace length=length+1 if sep==1
gen length1=length
expand length1, generate(newvar)
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
gen month=ysta+t

gen eqd=0
replace eqd=1 if month>=eq & eq!=0

gen unid=0
replace unid=1 if month>=unil & unil!=0

gen tit=0
replace tit=1 if month<eq 

gen com=0
replace com=1 if eq==0
replace com=0 if eq==0 & state=="Wisconsin" & month<(1986-1900)*12+1

*interaction terms
gen int1=unid*com
gen int2=unid*tit
gen int3=unid*eqd
*replace month=month+t
*****************************************************


*Some Additional variables
gen beg=round(COHAB1/12+1900)
gen birthh=round(A_3/12+1900)

egen yrr = cut(month), group(10)

bysort id:egen ini=sum(unid) 
sort id t
bysort id:replace ini=ini/t[_N]
gen see=0
replace see=1 if ini>0 & ini<1  
gen inter=see*unid
keep if see==0



*Basic
/*
reg faill unid inter see  i.tt i.beg i.A_4 i.birthh  if A_4<=56 , vce(cluster id)
mlogit fail2 unid inter see   t t2 t3 t4 i.beg i.A_4 i.birthh   if A_4<=56 , baseoutcome(0) difficult 
mlogit fail2 unid   i.tt i.beg i.A_4 i.birthh   if A_4<=56 & F_3==996, baseoutcome(0) difficult
mlogit fail2 change*unid   t t2 t3 t4  i.year i.A_4 i.birth EDUCAT i.RELIGION  if A_4<=56 , baseoutcome(0)
mlogit fail2 change*unid   t t2 t3 t4  i.year i.A_4 i.birth EDUCAT i.RELIGION  if A_4<=56 & (F_3==996 | (QDICMO-F_3)+36>month) , baseoutcome(0)
mlogit fail2 change*unid   t t2 t3 t4  i.year i.A_4 i.birth EDUCAT i.RELIGION  if A_4<=56 & F_3==996  , baseoutcome(0)

*Heterogeneity
mlogit fail2 int1 int2 int3 eqd com  t t2 t3 t4  i.year i.A_4 i.birth EDUCAT i.RELIGION  if A_4<=56 , baseoutcome(0)
mlogit fail2 int1 int2 int3 eqd com  t t2 t3 t4  i.year i.A_4 i.birth EDUCAT i.RELIGION  if A_4<=56 & (F_3==996 | (QDICMO-F_3)+36>month) , baseoutcome(0)
mlogit fail2 int1 int2 int3 eqd com  t t2 t3 t4  i.year i.A_4 i.birth EDUCAT i.RELIGION  if A_4<=56 & F_3==996 , baseoutcome(0)
*/


*Adapt for merge 
gen nsfh=0
gen nsfg=1
gen keep=0
replace keep=1 if F_3==996
rename birthh birth
gen female=1
rename AGE age
rename QDICMO itw
rename COHAB1 dati
gen agep=round(age-(itw-dati)/12)

*FINAL KEEP
keep fail2 unid inter see t beg state A_4 F_3 id age agep itw dati birth keep nsfh nsfg female int1 int2 int3 coll tit
keep if A_4<=56



append using "C:\Users\Fabio\Dropbox\JMP\empirical analysis\Analysis\duration\NSFHW1.dta"

*Rearranged stuff
egen tt= cut(t),group(20)
encode state,gen(stat)
egen beg1=cut(beg),group(5)
egen birth1=cut(birth),group(5)
ge t2=t^2
ge t3=t^3
ge t4=t^4
ge t5=t^5
gen censored=0
sort id t
bysort id:replace censored=1 if fail2[_N]==0
sort id t
gen sel=0
gen nid=_n
bysort id:replace sel=1 if _n==_N
1
bysort state:egen tot=mean(unid)

gen never=0
replace never=1 if tot==0
*drop if tot>=1 | tot<=0
drop if age<20 | age>=.

saveold "C:\Users\Fabio\Dropbox\JMP\empirical analysis\Analysis\duration\NSFG88.dta",replace


/*
********************
*Graphs Part
********************
bysort id: egen meanu=mean(unid)
*keep if  meanu==1 | meanu==0
sort id t
bysort id:egen max=max(t)
keep if max==t

*histograph of length
#delimit;
hist beg
         ,xlabel(1920(20)1988)
          xtitle(Year Cohabitation Began)
		  graphregion(color(white)) plotregion(color(white))
		   legend(row(1) label(1 "Simulated") ) ;
# delimit cr
graph export "C:\Users\Fabio\Dropbox\JMP\presentation\phd_apero_08_2019\hist_length.pdf", as(pdf) replace

*When relationship began from 1960 to 1985
keep if beg<=1985

gen time=.
replace time=1947.7 if beg>=1945 & beg<1955
replace time=1957.5 if beg>=1955 & beg<1960
replace time=1962.5 if beg>=1960 & beg<1965
replace time=1967.5 if beg>=1965 & beg<1970
replace time=1972.5 if beg>=1970 & beg<1975
replace time=1977.5 if beg>=1975 & beg<1980
replace time=1982.5 if beg>=1980 & beg<1985

/*
keep if beg<=1985 & beg>=1950
egen p0=min(beg)
egen p20=pctile(beg), p(20)
egen p40=pctile(beg), p(40)
egen p60=pctile(beg), p(60)
egen p80=pctile(beg), p(80)
egen p100=min(beg)

gen time=.
replace time=(p20+p0)/2 if beg<=p20
replace time=(p20+p40)/2 if beg<=p40 & beg>p20
replace time=(p40+p60)/2 if beg<=p60 & beg>p40
replace time=(p60+p80)/2 if beg<=p80 & beg>p60
replace time=(p80+p100)/2 if beg<=p100 & beg>p80
*/


bysort time never: egen tmean=mean(t)

#delimit;
twoway (line tmean time if never==0 & beg>=1965 & beg<=1985, color(gs8*.9) lwidth(0.45)) (line tmean time if never==1 &  beg>=1965 & beg<=1985,color(maroon*.75) lwidth(0.45)),
         graphregion(color(white)) plotregion(color(white)) xtitle(Cohabitation Began - Year) ytitle(Cohabitation length (Months))
         legend(order(1 2) row(1) label(1 "Became Unilateral States") label(2 "Never Unilateral States") 
	label(3 "Data") );
# delimit cr
graph export "C:\Users\Fabio\Dropbox\JMP\presentation\phd_apero_08_2019\ev_length.pdf", as(pdf) replace
/*
mlogit fail2 unid   i.tt i.beg i.stat i.birth  , baseoutcome(0) difficult
mlogit fail2 unid   i.tt i.beg i.stat i.birth    if keep==1, baseoutcome(0) difficult
mlogit fail2 unid   t i.beg1 i.stat i.birth1  , baseoutcome(1) difficult

gen censored=1
replace censored=0 if fail2==1
replace censored=0 if fail2==2
bysort id:gen sel=1 if _n==_N
keep if sel==1
tabulate stat, generate(s)
ge birth2=birth^2
gen birth3=birth^3
tabulate beg, generate(be)
cnreg t unid beg s01-s24 b1-b27 if sel==1, censored(censored)


*rearrange the data
gen choice=3
expand(choice)

*gen choice variable

gen y=0
bysort id t:ge cho=_n

replace y=1 if cho==1 & fail2==0
replace y=1 if cho==2 & fail2==1
replace y=1 if cho==3 & fail2==2

ge rel2=0
replace rel2=1 if cho==3 


nlogitgen rel = cho(3, 1 2)
egen nid=group(id t)
nlogittree cho rel, choice(y) case(nid)


cmset id t y

 nlogit y  i.beg i.stat i.birth || rel:i.tt || cho:,case(nid)         
