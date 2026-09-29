*! contdid_simdata.do — generate simulated example data
*! DGP: two-period panel, continuous dose D, true ATT(d) = 0.5*d
*! Run from the contdid/ package root.

version 16
clear all
set seed 12345

cd "/Users/niuhaoyu/Documents/open code/contdid"

set obs 4000
gen id = _n
* ~20% untreated (D=0), ~80% treated with dose ~ U(0,1)
gen d = cond(runiform() > 0.2, runiform(), 0)

gen mu = 2 + rnormal(0, 1)      // unit fixed effect
gen e0 = rnormal(0, 0.3)
gen e1 = rnormal(0, 0.3)
gen y0 = mu + e0                 // pre-treatment outcome
gen y1 = mu + 0.5*d + e1         // post-treatment outcome; true ATT(d) = 0.5*d

gen t = 0
gen y = y0
tempfile pre
save `pre'

replace t = 1
replace y = y1
append using `pre'

keep id t y d
order id t y d
sort id t

save "data/contdid_sim.dta", replace
