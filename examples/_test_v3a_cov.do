*! _test_v3a_cov.do - v3a covariates test (conditional parallel trends)
*! True ATT(d) = (0.5 + 0.2 * mean(x|D>0)) * d
version 16
clear all
set more off
adopath + "."

use "data/contdid_cov_sim.dta", clear

qui summarize x if d > 0
local xbar = r(mean)
di as result "xbar(treated) = " %8.5f `xbar'

* (1) unconditional estimator (expected BIASED)
contdid y, unit(id) time(t) dose(d) degree(1) nknots(0) npoints(9) reps(0)
matrix u = r(attd)

* (2) covariates estimator (expected ~unbiased)
contdid y, unit(id) time(t) dose(d) covariates(x) degree(1) nknots(0) npoints(9) reps(0)
matrix c = r(attd)

di as text _n "   d       true      uncond    |bias|     cond      |bias|"
local maxbu = 0
local maxbc = 0
forvalues k = 1/9 {
    local dk  = c[`k',1]
    local tru = (0.5 + 0.2*`xbar') * `dk'
    local bu  = abs(u[`k',2] - `tru')
    local bc  = abs(c[`k',2] - `tru')
    if `bu' > `maxbu' local maxbu = `bu'
    if `bc' > `maxbc' local maxbc = `bc'
    di as text %6.4f `dk' "  " %9.6f `tru' "  " %9.6f u[`k',2] "  " %8.5f `bu' "  " %9.6f c[`k',2] "  " %8.5f `bc'
}
di as result _n "max|bias| unconditional = " %8.5f `maxbu'
di as result    "max|bias| covariates    = " %8.5f `maxbc'
assert `maxbc' < `maxbu'
assert `maxbc' < 0.05
di as result "V3A TEST PASS"
