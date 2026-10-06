*! v3a_simdata.do - conditional parallel trends DGP for contdid v3a
*! True ATT(d) = (0.5 + 0.2 * mean(x|D>0)) * d ; unconditional PT FAILS.
*! Run from the contdid/ package root.

version 16
clear all
set seed 20261006

set obs 3000
gen long id = _n
gen double x = rnormal(0, 1)

* treatment selection depends on x  ->  E[x | D=d] varies with d
gen double p = invlogit(0.8 * x)
gen byte treated = runiform() < p
gen double d = 0
replace d = runiform() if treated          // dose ~ U(0,1) among treated

* untreated trend depends on x  ->  unconditional parallel trends fails
* treated outcome change: (0.5 + 0.2*x) * d  ->  ATT_x(d) = (0.5+0.2x) d
gen double mu  = rnormal(0, 0.5)
gen double e0  = rnormal(0, 0.3)
gen double e1  = rnormal(0, 0.3)
gen double dy  = 0.3 * x + (0.5 + 0.2 * x) * d + e1 - e0
gen double y0  = mu + e0
gen double y1  = mu + dy

gen byte t = 0
gen double y = y0
tempfile pre
save `pre'

replace t = 1
replace y = y1
append using `pre'

keep id t y d x
order id t y d x
sort id t

save "data/contdid_cov_sim.dta", replace
export delimited id t y d x using "examples/reference/contdid_cov_sim.csv", replace
display "wrote data/contdid_cov_sim.dta and examples/reference/contdid_cov_sim.csv"
