*! _test_coverage.do — cluster bootstrap CI should have ≈95% coverage
version 16
clear all
cd "/Users/niuhaoyu/Documents/open code/contdid"
adopath + "/Users/niuhaoyu/Documents/open code/contdid"

local niter   = 100
local covered = 0

forvalues it = 1/`niter' {
    * small-sample DGP, true ATT(d) = 0.5*d
    clear
    set seed `it'
    set obs 200
    gen id = _n
    gen d  = cond(runiform() > 0.2, runiform(), 0)
    gen mu = 2 + rnormal(0, 1)
    gen e0 = rnormal(0, 0.3)
    gen e1 = rnormal(0, 0.3)
    gen y0 = mu + e0
    gen y1 = mu + 0.5*d + e1
    gen t  = 0
    gen y  = y0
    tempfile pre
    save `pre'
    replace t = 1
    replace y = y1
    append using `pre'
    keep id t y d

    qui contdid y, unit(id) time(t) dose(d) npoints(11) reps(199) seed(`it')
    matrix b = r(attd)
    * point 6 is the midpoint d≈0.5, true ATT≈0.25
    local lb = b[6,4]
    local ub = b[6,5]
    if `lb' <= 0.25 & 0.25 <= `ub' {
        local covered = `covered' + 1
    }
}

local cov = `covered' / `niter'
di as text "COVERAGE = " as result `cov' "  (target ≈ 0.95)"
assert `cov' > 0.85 & `cov' < 1.0
di as result "PASS: coverage in expected range"
