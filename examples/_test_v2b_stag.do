*! staggered via contdid command: true ATT(d)=0.5d
version 16
clear all
cd "/Users/niuhaoyu/Documents/open code/contdid"
adopath + "/Users/niuhaoyu/Documents/open code/contdid/ado"

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

contdid y, unit(id) time(t) dose(d) gvar(g) npoints(5) reps(0)

matrix b = r(attd)
matrix a = r(acrt)
local maxerr = 0
forvalues i = 1/5 {
    local dk = b[`i',1]
    local truth = 0.5 * `dk'
    local err = abs(b[`i',2] - `truth')
    if `err' > `maxerr' local maxerr = `err'
}
di as text "STAGGERED ATT max err = " as result `maxerr'
assert `maxerr' < 0.05
di as result "PASS: staggered via contdid recovers ATT(d)=0.5d"
