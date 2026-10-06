*! _test_v3a_cband_cov.do - joint coverage of the covariate uniform band
version 16
clear all
set more off
adopath + "."

local R = 100
local cover = 0
forvalues r = 1/`R' {
    qui {
        clear
        set seed `= 10000 + `r''
        set obs 1200
        gen long id = _n
        gen double x = rnormal(0, 1)
        gen double p = invlogit(0.8 * x)
        gen byte treated = runiform() < p
        gen double d = 0
        replace d = runiform() if treated
        gen double mu = rnormal(0, 0.5)
        gen double e0 = rnormal(0, 0.3)
        gen double e1 = rnormal(0, 0.3)
        gen double dy = 0.3 * x + (0.5 + 0.2 * x) * d + e1 - e0
        gen double y0 = mu + e0
        gen double y1 = mu + dy
        gen byte t = 0
        gen double y = y0
        tempfile pre
        save `pre'
        replace t = 1
        replace y = y1
        append using `pre'
        summarize x if d > 0
        local xbar = r(mean)
        contdid y, unit(id) time(t) dose(d) covariates(x) degree(1) nknots(0) npoints(9) cband reps(200) seed(1)
        matrix cc = r(cb_att)
        local ok = 1
        forvalues k = 1/9 {
            local dk  = cc[`k',1]
            local tru = (0.5 + 0.2 * `xbar') * `dk'
            if `tru' < cc[`k',3] | `tru' > cc[`k',4] local ok = 0
        }
    }
    if `ok' local cover = `cover' + 1
}
di as result _n "JOINT COVERAGE (nominal 0.95) = " `cover'/`R'
assert `cover'/`R' >= 0.85
di as result "V3A CBAND COVERAGE TEST PASS"
