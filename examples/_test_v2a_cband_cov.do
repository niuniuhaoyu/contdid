*! uniform-band coverage Monte Carlo (target ~95%)
version 16
clear all
cd "/Users/niuhaoyu/Documents/open code/contdid"
adopath + "/Users/niuhaoyu/Documents/open code/contdid"

local niter   = 50
local covered = 0

forvalues it = 1/`niter' {
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

    qui contdid y, unit(id) time(t) dose(d) npoints(11) reps(499) seed(`it') cband
    matrix cb = r(cb_att)
    local allcov = 1
    forvalues k = 1/11 {
        local dk = cb[`k',1]
        local lb = cb[`k',3]
        local ub = cb[`k',4]
        if `lb' > 0.5*`dk' | `ub' < 0.5*`dk' local allcov = 0
    }
    if `allcov' == 1 local covered = `covered' + 1
}

local cov = `covered' / `niter'
di as text "UNIFORM BAND COVERAGE = " as result `cov' "  (target ~0.95)"
assert `cov' > 0.85 & `cov' < 1.0
di as result "PASS: uniform band coverage in expected range"
