* export_attd.do - export contdid r(attd) for cross-check vs base-R / R contdid
version 16
clear all
set more off
adopath + "."

use "data/contdid_sim.dta", clear

contdid y, unit(id) time(t) dose(d) degree(1) nknots(0) npoints(9) reps(50) seed(12345)

display "=== return list ==="
return list

matrix b = r(attd)
display "r(attd) rows=" rowsof(b) " cols=" colsof(b)
mat list b

clear
svmat double b
export delimited using "examples/reference/stata_attd.csv", replace
display "wrote examples/reference/stata_attd.csv"
