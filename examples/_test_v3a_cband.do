*! _test_v3a_cband.do - uniform band with covariates
version 16
clear all
set more off
adopath + "."

use "data/contdid_cov_sim.dta", clear
contdid y, unit(id) time(t) dose(d) covariates(x) degree(1) nknots(0) npoints(9) cband reps(200) seed(12345)

matrix c  = r(attd)
matrix cc = r(cb_att)
assert rowsof(cc) == 9
forvalues k = 1/9 {
    assert abs(cc[`k',2] - c[`k',2]) < 1e-8
    assert cc[`k',3] <= c[`k',2] + 1e-9
    assert cc[`k',4] >= c[`k',2] - 1e-9
    assert cc[`k',4] > cc[`k',3]
}
di as result "crit_att = " %7.4f r(crit_att)
di as result "V3A CBAND TEST PASS"
