*! _test_recover.do — verify contdid recovers true ATT(d) = 0.5*d
version 16
clear all
cd "/Users/niuhaoyu/Documents/open code/contdid"
adopath + "/Users/niuhaoyu/Documents/open code/contdid/ado"
use "data/contdid_sim.dta", clear

contdid y, unit(id) time(t) dose(d) npoints(5) seed(1)
matrix b = r(attd)

local maxerr = 0
forvalues i = 1/5 {
    local dk = b[`i',1]
    local att = b[`i',2]
    local err = abs(`att' - 0.5*`dk')
    if `err' > `maxerr' local maxerr = `err'
}
di as text "MAX ERR = " as result `maxerr'
assert `maxerr' < 0.05
di as result "PASS: contdid recovers true ATT(d)"
