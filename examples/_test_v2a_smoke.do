*! smoke test: v2a runs on default (degree=1 linear), recover 0.5d
version 16
clear all
cd "/Users/niuhaoyu/Documents/open code/contdid"
adopath + "/Users/niuhaoyu/Documents/open code/contdid/ado"
use "data/contdid_sim.dta", clear

contdid y, unit(id) time(t) dose(d) npoints(5) reps(0)
matrix b = r(attd)
matrix a = r(acrt)
local maxerr = 0
forvalues i = 1/5 {
    local dk = b[`i',1]
    local att = b[`i',2]
    local err = abs(`att' - 0.5*`dk')
    if `err' > `maxerr' local maxerr = `err'
}
di as text "ATT MAX ERR (degree=1) = " as result `maxerr'
assert `maxerr' < 0.05
* ACRT should be ~ constant 0.5 (linear slope)
local acrt1 = a[1,2]
local acrt5 = a[5,2]
di as text "ACRT range = " as result `acrt1' " .. " `acrt5'
assert abs(`acrt1' - 0.5) < 0.05 & abs(`acrt5' - 0.5) < 0.05
di as result "PASS: v2a degree(1) recovers 0.5d and ACRT=0.5"
