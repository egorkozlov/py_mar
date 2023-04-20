***************************************************************
*Get people of the 45-55 birth cohort cohabiting in this wave
***************************************************************
clear all

*Import data
use "C:\Users\Fabio\Dropbox\JMP\empirical analysis\NSFH\Wave 2\ICPSR_06906\DS0001\06906-0001-Data.dta" 

*Keep of the right birth cohorts
keep if MA4>=45 & MA4<55 

*Keep if cohabiting in the wave
*keep if MSCALCA==1

*Keep age only as a variable
keep MA8

rename MA8 age

*Save
export delimited using "C:\Users\Fabio\Dropbox\py_mar1\py_mar\age_drop.csv",replace
