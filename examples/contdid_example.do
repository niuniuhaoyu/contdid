*! contdid_example.do — one-click reproducible example
*! Run from the contdid/ package root.

version 16
clear all
cd "/Users/niuhaoyu/Documents/open code/contdid"
adopath + "/Users/niuhaoyu/Documents/open code/contdid/ado"

* 1. Load the simulated data (true dose-response ATT(d) = 0.5*d)
*    Regenerate with: do "examples/contdid_simdata.do"
use "data/contdid_sim.dta", clear

* 2. Estimate the dose-response function
contdid y, unit(id) time(t) dose(d) npoints(20) reps(999) seed(12345) graph

* 3. Save the figure
graph export "examples/contdid_dose_response.png", width(1200) replace

di as result _n "Done. See examples/contdid_dose_response.png"
