*! nonlinear recovery: quadratic DGP, degree(2), ATT(d)=0.3d+0.4d^2, ACRT=0.3+0.8d
version 16
clear all
cd "/Users/niuhaoyu/Documents/open code/contdid"
adopath + "/Users/niuhaoyu/Documents/open code/contdid"

* quadratic DGP
clear
set seed 7
set obs 4000
gen id = _n
gen d = cond(runiform() > 0.2, runiform(), 0)
gen mu = 2 + rnormal(0,1)
gen e0 = rnormal(0,0.2)
gen e1 = rnormal(0,0.2)
gen y0 = mu + e0
gen y1 = mu + 0.3*d + 0.4*d^2 + e1
gen t = 0
gen y = y0
tempfile pre
save `pre'
replace t = 1
replace y = y1
append using `pre'
keep id t y d

contdid y, unit(id) time(t) dose(d) npoints(11) reps(199) seed(1) degree(2)

matrix b = r(attd)
matrix a = r(acrt)
local maxerr_att = 0
local maxerr_acr = 0
forvalues i = 1/11 {
    local dk = b[`i',1]
    local att = b[`i',2]
    local truth_att = 0.3*`dk' + 0.4*`dk'^2
    local err1 = abs(`att' - `truth_att')
    if `err1' > `maxerr_att' local maxerr_att = `err1'
    local acr = a[`i',2]
    local truth_acr = 0.3 + 0.8*`dk'
    local err2 = abs(`acr' - `truth_acr')
    if `err2' > `maxerr_acr' local maxerr_acr = `err2'
}
di as text "ATT max err (quadratic) = " as result `maxerr_att'
di as text "ACRT max err (quadratic) = " as result `maxerr_acr'
assert `maxerr_att' < 0.05
assert `maxerr_acr' < 0.15
* CI sanity: lb <= att <= ub at midpoint (i=6)
assert b[6,4] <= b[6,2] & b[6,2] <= b[6,5]
di as result "PASS: degree(2) recovers quadratic ATT(d) and ACRT(d)"
