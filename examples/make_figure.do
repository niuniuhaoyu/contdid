*! generate a polished figure for sharing: nonlinear dose-response ATT(d) + ACRT(d)
version 16
clear all
cd "/Users/niuhaoyu/Documents/open code/contdid"
adopath + "/Users/niuhaoyu/Documents/open code/contdid/ado"
set scheme s1color

* ---- DGP: concave dose-response with diminishing returns ----
clear
set seed 2026
set obs 3000
gen id = _n
gen d  = cond(runiform() > 0.2, runiform(), 0)     // ~20% untreated
gen mu = 2 + rnormal(0, 1)
gen e0 = rnormal(0, 0.25)
gen e1 = rnormal(0, 0.25)
gen y0 = mu + e0
gen y1 = mu + (0.8*d - 0.3*d^2) + e1               // true ATT(d) = 0.8d - 0.3d^2
gen t  = 0
gen y  = y0
tempfile pre
save `pre'
replace t = 1
replace y = y1
append using `pre'
keep id t y d

* ---- estimate: cubic B-spline, uniform confidence band ----
contdid y, unit(id) time(t) dose(d) degree(3) nknots(2) npoints(30) reps(299) seed(1) cband

* ---- build plot data ----
preserve
clear
matrix a = r(attd)
matrix cb = r(cb_att)
matrix ac = r(acrt)
local np = rowsof(a)
set obs `np'
gen double dd  = .
gen double att = .
gen double lb  = .
gen double ub  = .
gen double cbl = .
gen double cbu = .
gen double acrt = .
forvalues i = 1/`np' {
    qui replace dd   = a[`i',1]  in `i'
    qui replace att  = a[`i',2]  in `i'
    qui replace lb   = a[`i',4]  in `i'
    qui replace ub   = a[`i',5]  in `i'
    qui replace cbl  = cb[`i',3] in `i'
    qui replace cbu  = cb[`i',4] in `i'
    qui replace acrt = ac[`i',2] in `i'
}

* ---- figure: ATT(d) with uniform band (top) and ACRT(d) (bottom) ----
twoway (rarea cbu cbl dd, color(navy%12)) ///
       (line att dd, lcolor(navy) lwidth(medthick)), ///
    yline(0, lcolor(gs8) lpattern(dash)) ///
    legend(off) ///
    xtitle("Treatment dose (d)", size(medlarge)) ///
    ytitle("ATT(d)", size(medlarge)) ///
    title("Dose-response: ATT(d)", size(large)) ///
    subtitle("with 95% uniform confidence band", size(medsmall) color(gs6)) ///
    note("B-spline (degree 3, 2 knots) · contdid · Callaway–Goodman-Bacon–Sant'Anna (2024)", size(small) color(gs6)) ///
    name(att_plot, replace)

twoway (line acrt dd, lcolor(maroon) lwidth(medthick)) ///
    (rarea cbu cbl dd, color(navy%6)), ///
    yline(0, lcolor(gs8) lpattern(dash)) ///
    legend(off) ///
    xtitle("Treatment dose (d)", size(medlarge)) ///
    ytitle("ACRT(d)", size(medlarge)) ///
    title("Average causal response: ACRT(d)", size(large)) ///
    subtitle("marginal effect of an additional unit of dose", size(medsmall) color(gs6)) ///
    name(acrt_plot, replace)

graph combine att_plot acrt_plot, cols(1) xsize(8) ysize(9) ///
    title("Continuous-treatment difference-in-differences", size(large)) ///
    name(combined, replace)

graph export "examples/contdid_figure.png", width(2000) replace

restore
di as result _n "Figure exported: examples/contdid_figure.png"
