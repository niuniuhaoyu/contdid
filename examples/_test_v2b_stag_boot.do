*! staggered with bootstrap: CI non-degenerate, covers truth
version 16
clear all
cd "/Users/niuhaoyu/Documents/open code/contdid"
adopath + "/Users/niuhaoyu/Documents/open code/contdid"

clear
set seed 42
set obs 3000
gen id = _n
gen d = cond(runiform() > 0.3, runiform(), 0)
gen g = 0
replace g = 2 + floor(3*runiform()) if d > 0
gen mu = 1 + rnormal(0,1)
forvalues t = 1/4 {
    gen y`t' = mu + 0.5*d*(`t' >= g) + rnormal(0,0.3)
}
reshape long y, i(id) j(t)

contdid y, unit(id) time(t) dose(d) gvar(g) npoints(5) reps(99) seed(1)

matrix b = r(attd)
* CI sanity: at each point, lb <= att <= ub and se > 0
forvalues i = 1/5 {
    assert b[`i',3] > 0
    assert b[`i',4] <= b[`i',2] & b[`i',2] <= b[`i',5]
}
* midpoint d~0.5 truth~0.25 covered by CI
assert b[3,4] <= 0.25 & 0.25 <= b[3,5]
di as result "PASS: staggered bootstrap CI valid, midpoint covers truth"
