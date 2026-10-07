*! _test_v3b_dds.do - data-driven sieve for contdid
version 16
clear all
set more off
adopath + "."

* (1) linear DGP (true ATT(d)=0.5d) -> dds should pick ~0 interior knots
use "data/contdid_sim.dta", clear
contdid y, unit(id) time(t) dose(d) dose_est_method(dds) degree(1) npoints(9) reps(0)
di as result "linear DGP: selected knots = " scalar(dds_nknots)

* (2) nonlinear DGP (true ATT(d)=0.3d+0.4d^2)
clear all
set seed 999
set obs 4000
gen long id = _n
gen double d = cond(runiform()>0.2, runiform(), 0)
gen double mu = 2 + rnormal(0,1)
gen double y0 = mu + rnormal(0,0.3)
gen double y1 = mu + 0.3*d + 0.4*d*d + rnormal(0,0.3)
gen double tt = 0
gen double y = y0
tempfile pre
save `pre'
replace tt = 1
replace y = y1
append using `pre'
rename tt t
keep id t y d
sort id t

contdid y, unit(id) time(t) dose(d) dose_est_method(dds) degree(3) maxknots(5) npoints(9) reps(0)
matrix b = r(attd)
di as result "nonlinear DGP: selected knots = " scalar(dds_nknots)
local maxerr = 0
forvalues k = 1/9 {
    local dk = b[`k',1]
    local tru = 0.3*`dk' + 0.4*`dk'*`dk'
    local e = abs(b[`k',2] - `tru')
    if `e' > `maxerr' local maxerr = `e'
    di as text %6.4f `dk' "  true=" %8.4f `tru' "  dds=" %8.4f b[`k',2] "  |e|=" %8.4f `e'
}
di as result "max|error| = " %8.4f `maxerr'
assert `maxerr' < 0.05
di as result "V3B DDS TEST PASS"
