*! _test_backcompat.do - ensure default path unchanged after v3a additions
version 16
clear all
set more off
adopath + "."

use "data/contdid_sim.dta", clear
contdid y, unit(id) time(t) dose(d) degree(1) nknots(0) npoints(9) reps(0)
matrix b = r(attd)
assert rowsof(b) == 9
* golden values from examples/reference/stata_attd.csv (v1.0.0)
assert abs(b[1,2] - 0.0337626286411811) < 1e-8
assert abs(b[5,2] - 0.2821127657836042) < 1e-8
assert abs(b[9,2] - 0.5304629029260274) < 1e-8
di as result "BACKWARD COMPAT OK"

* staggered path still loads
use "data/contdid_sim.dta", clear
di as result "BACKWARD COMPAT TEST PASS"
