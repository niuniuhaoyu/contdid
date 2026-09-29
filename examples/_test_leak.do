*! quick leak check: call contdid twice
version 16
clear all
cd "/Users/niuhaoyu/Documents/open code/contdid"
adopath + "/Users/niuhaoyu/Documents/open code/contdid"
use "data/contdid_sim.dta", clear

contdid y, unit(id) time(t) dose(d) npoints(5) reps(99) seed(1)
di as result "--- second call ---"
contdid y, unit(id) time(t) dose(d) npoints(5) reps(99) seed(2)
di as result "PASS: no preserve leak"
