*! contdid_example.do — one-click reproducible example
*! Run from the contdid/ package root.

version 16
clear all

* 1. Load the simulated data (true dose-response ATT(d) = 0.5*d)
*    Regenerate with: do "examples/contdid_simdata.do"
use "data/contdid_sim.dta", clear

* 2. Make the command available (skip if already installed)
capture which contdid
if _rc {
    adopath + "."
}

* 3. Estimate the dose-response (cubic spline + uniform band + plot)
contdid y, unit(id) time(t) dose(d) degree(3) nknots(2) npoints(20) reps(999) seed(12345) cband graph

* 4. Save the figure
graph export "examples/contdid_dose_response.png", width(1200) replace

di as result _n "Done. See examples/contdid_dose_response.png"
