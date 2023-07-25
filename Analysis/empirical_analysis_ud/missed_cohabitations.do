*****************************************************
**Build first Cohabitation Dataset
*****************************************************

*****************************************************
**1 NSFH
*****************************************************
clear  all
use "C:\Users\Fabio\Dropbox\JMP\empirical analysis\NSFH\Wave 1\ICPSR_06041\DS0001\06041-0001-Data.dta" 


*Observations with inconsistencies
gen inconsistency=0

*How many cohabitation I am missing?
gen cohabitation_num=0

**************************************************************************************
*First consider those currently cohabiting and never married ( BKMK2==1 &  M2CP01==5)
replace cohabitation_num=1+M145 if M145<=6 & BKMK2==1 &  M2CP01==5


*The next line OR the other following two lines should be commented
*replace inconsistency=1 if M145>6 & BKMK2==1 &  M2CP01==5
*replace cohabitation_num=1 if M145>6 & BKMK2==1 &  M2CP01==5 & M140<=1 
*replace inconsistency=1 if M145>6 & BKMK2==1 &  M2CP01==5 & M140>1
***************************************************************************************


**************************************************************************************
*Consider those not currently cohabiting and never married ( BKMK2==2 &  M2CP01==5)

replace cohabitation_num=M140 if M140<=7 & BKMK2==2 &  M2CP01==5
*replace inconsistency=1 if M145>7 & BKMK2==2 &  M2CP01==5
***************************************************************************************


**************************************************************************************
*Cohabitation with first spouse? To ever married M2CP01!=5 
************************************************************************************
replace cohabitation_num=1 if M113==1 & M2CP01!=5
*replace inconsistency=1 if M2CP01==5 & M113<6


**************************************************************************************
*Cohabitations before first spouse To ever married M2CP01!=5 
************************************************************************************
replace cohabitation_num=cohabitation_num+M117 if  M2CP01!=5 & M117<=7 





**************************************************************************************
*Did you cohabit with second wife?. People married at least twice
* M2CP01!=5 & M103NUM>=1 & M103NUM<=4
************************************************************************************
replace cohabitation_num=cohabitation_num+1 if  M2CP01!=5 &  M103NUM>=1 & M103NUM<=4 & M122==1


**************************************************************************************
*Cohabitations between first and second marriage . People married at least twice
* M2CP01!=5 & M103NUM>=1 & M103NUM<=4
************************************************************************************
replace cohabitation_num=cohabitation_num+M124 if  M2CP01!=5 &  M103NUM>=1 & M103NUM<=4 & M124<=6


**************************************************************************************
*Cohabitations after second marriage if no more married and not currently cohabiting (M95==2) & (BKMK2==2)
************************************************************************************
replace cohabitation_num=cohabitation_num+M131 if  M95<=2 & M95>0 & BKMK2==2 & M131<=7


**************************************************************************************
*Cohabitations after second marriage if no more married and currently cohabiting (M95==1) & (BKMK2==2)
************************************************************************************
replace cohabitation_num=cohabitation_num+1 if  M95<=2  & M95>0 & BKMK2==1 
replace cohabitation_num=cohabitation_num+M136 if  M95<=2  & M95>0 & BKMK2==1 & M136<=5

gen diff=cohabitation_num- NUMCOHAB

*Now see for how many observations we lost 1 (2) cohabitations
gen pdiff=0
replace pdiff=1 if diff>=1
sum pdiff if M95<=2 [weight=SAMWT]

gen pdiff2=0
replace pdiff2=1 if diff>=2
sum pdiff2 if M95<=2 [weight=SAMWT]
