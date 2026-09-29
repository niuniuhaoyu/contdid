*! cband smoke: runs and produces crit values + bands
version 16
clear all
cd "/Users/niuhaoyu/Documents/open code/contdid"
adopath + "/Users/niuhaoyu/Documents/open code/contdid"
use "data/contdid_sim.dta", clear

contdid y, unit(id) time(t) dose(d) npoints(11) reps(999) seed(1) cband
di as text "crit_att = " as result r(crit_att) as text ", crit_acrt = " as result r(crit_acrt)
* uniform band should be wider than pointwise (crit > 1.96 approx)
assert r(crit_att) > 1.5 & r(crit_att) < 6
assert r(crit_acrt) > 1.5 & r(crit_acrt) < 6
* band must cover the point estimate
matrix cb = r(cb_att)
assert cb[6,3] <= cb[6,2] & cb[6,2] <= cb[6,4]
di as result "PASS: cband runs, crit ~ " r(crit_att) " / " r(crit_acrt)
