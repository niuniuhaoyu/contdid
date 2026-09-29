*! verify portable example works from package root
version 16
clear all
cd "/Users/niuhaoyu/Documents/open code/contdid"
use "data/contdid_sim.dta", clear
capture which contdid
if _rc {
    adopath + "ado"
}
contdid y, unit(id) time(t) dose(d) npoints(5) reps(99) seed(1)
di as result "PASS: portable paths resolve from package root"
