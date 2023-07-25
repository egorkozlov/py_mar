*****************************************
*File to run all file linked to all the
*empirical results linked to unilateral
*divorce and partnership choice
*Wired commands: savesome
*****************************************
clear all
cd "C:\Users\Fabio\Dropbox\py_mar_3\Analysis\empirical_analysis_ud"

capture log close
log using "C:\Users\Fabio\Dropbox\py_mar_3\Analysis\unid_all", text replace
capture drop _all


*Cohabitation/marriage choice - NSFH
do relationship_nsfh_new.do
do relationship_nsfh_smc.do

*Cohabitation/marriage choice - NSFG
*do relationship_nsfg_new.do

*Cohabitation duration - NSFH
do duration_nsfh_new.do

*Cohabitation duration - NSFG
*do duration_nsfg.do

log close