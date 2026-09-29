*! export data for R cross-check
version 16
clear all
cd "/Users/niuhaoyu/Documents/open code/contdid"
use "data/contdid_sim.dta", clear
export delimited using "examples/reference/contdid_sim.csv", replace
di "exported"
