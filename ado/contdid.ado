*! contdid: Difference-in-Differences with a Continuous Treatment
*! version 0.3.0  2026-09-29  Haoyu Niu
*! Implements Callaway, Goodman-Bacon & Sant'Anna (2024/2025)
*! Two-period or staggered continuous-treatment DiD. ATT(d) and ACRT(d) via B-spline (default linear).

program define contdid, rclass
    version 16

    syntax varlist(max=1 numeric) [if] [in], ///
        unit(varname numeric) ///           individual id
        time(varname numeric) ///           time period (two periods)
        dose(varname numeric) ///           continuous treatment dose
        [npoints(integer 20) ///             number of evaluation doses
         level(real 95) ///                  confidence level (%)
         seed(integer 12345) ///             RNG seed (bootstrap)
         reps(integer 999) ///               bootstrap replications
         cluster(varname) ///                cluster variable (default = unit)
         degree(integer 1) ///               B-spline degree (1=linear)
         knots(numlist) ///                  explicit interior knots
         nknots(integer 0) ///               number of interior knots (quantile)
         gvar(varname numeric) ///           treatment timing (0=never treated; enables staggered)
         cband ///                            uniform confidence band (sup-t)
         GRaph]                               // dose-response plot

    local depvar `varlist'

    marksample touse
    preserve
    qui keep if `touse'

    * ---------- validation ----------
    qui summarize `dose', meanonly
    if r(min) < 0 {
        di as error "contdid: dose variable `dose' must be non-negative"
        exit 198
    }
    qui levelsof `time', local(tvals)
    local ntvals : word count `tvals'
    if "`gvar'" == "" & `ntvals' != 2 {
        di as error "contdid: supports only two time periods; found `ntvals' distinct values"
        exit 198
    }
    qui count if `dose' == 0
    if r(N) == 0 {
        di as error "contdid: no untreated units (dose==0); ATT(d) baseline requires an untreated group"
        exit 198
    }
    if `npoints' < 2 {
        di as error "contdid: npoints() must be at least 2"
        exit 198
    }
    if `degree' < 1 {
        di as error "contdid: degree() must be >= 1"
        exit 198
    }
    if "`knots'" != "" & `nknots' > 0 {
        di as error "contdid: specify either knots() or nknots(), not both"
        exit 198
    }
    if "`cband'" != "" & `reps' == 0 {
        di as error "contdid: cband requires reps() > 0 (multiplier bootstrap draws)"
        exit 198
    }

    * ---------- load Mata functions (once per session) ----------
    capture mata: bspline_basis(J(1,1,0.5), 1, J(0,1,0))
    if _rc {
        findfile "bspline.mata"
        qui do "`r(fn)'"
    }

    * ---------- staggered adoption path ----------
    if "`gvar'" != "" {
        qui summarize `gvar', meanonly
        if r(min) < 0 {
            di as error "contdid: gvar() must be >= 0 (0 = never treated)"
            exit 198
        }
        qui levelsof `time', local(_tvals)
        local _T : word count `_tvals'
        matrix _tvalsm = J(1, `_T', 0)
        local _ti = 0
        foreach _tv of local _tvals {
            local _ti = `_ti' + 1
            matrix _tvalsm[1, `_ti'] = `_tv'
        }
        qui summarize `dose' if `dose' > 0, meanonly
        local dmin = r(min)
        local dmax = r(max)
        local knotmatname ""
        if "`knots'" != "" {
            matrix _knots = (`knots')
            local knotmatname "_knots"
        }
        qui reshape wide `depvar', i(`unit') j(`time')
        local _yvars ""
        foreach _tv of local _tvals {
            local _yvars "`_yvars' `depvar'`_tv'"
        }
        tempname attm acrm
        matrix `attm' = J(`npoints', 5, .)
        matrix `acrm' = J(`npoints', 5, .)
        mata: _contdid_stag_run("`dose'", "`gvar'", "`_yvars'", "_tvalsm", `degree', "`knotmatname'", `nknots', `npoints', `dmin', `dmax', "_att", "_acrt")
        forvalues k = 1/`npoints' {
            local dk = `dmin' + (`dmax' - `dmin') * (`k' - 1) / (`npoints' - 1)
            matrix `attm'[`k',1] = `dk'
            matrix `attm'[`k',2] = _att[1,`k']
            matrix `acrm'[`k',1] = `dk'
            matrix `acrm'[`k',2] = _acrt[1,`k']
        }

        * staggered cluster bootstrap
        if `reps' > 0 {
            set seed `seed'
            if "`cluster'" != "" {
                local bscmd "bsample, cluster(`cluster')"
            }
            else {
                local bscmd "bsample"
            }
            tempfile est
            qui save `est'
            matrix boot = J(`reps', `npoints', .)
            matrix boota = J(`reps', `npoints', .)
            forvalues b = 1/`reps' {
                qui use `est', clear
                qui `bscmd'
                mata: _contdid_stag_run("`dose'", "`gvar'", "`_yvars'", "_tvalsm", `degree', "`knotmatname'", `nknots', `npoints', `dmin', `dmax', "_batt", "_bacrt")
                forvalues k = 1/`npoints' {
                    matrix boot[`b',`k']  = _batt[1,`k']
                    matrix boota[`b',`k'] = _bacrt[1,`k']
                }
            }
            qui use `est', clear
            local plo = (100 - `level') / 2
            local phi = 100 - `plo'
            svmat boot, names(b_)
            forvalues k = 1/`npoints' {
                qui summarize b_`k'
                matrix `attm'[`k',3] = r(sd)
                qui centile b_`k', centile(`plo' `phi')
                matrix `attm'[`k',4] = r(c_1)
                matrix `attm'[`k',5] = r(c_2)
            }
            svmat boota, names(c_)
            forvalues k = 1/`npoints' {
                qui summarize c_`k'
                matrix `acrm'[`k',3] = r(sd)
                qui centile c_`k', centile(`plo' `phi')
                matrix `acrm'[`k',4] = r(c_1)
                matrix `acrm'[`k',5] = r(c_2)
            }
        }

        matrix colnames `attm' = d ATT se lb ub
        matrix colnames `acrm' = d ACRT se lb ub
        di as text _n "Dose-response ATT(d)  (staggered; B-spline degree `degree'; " ///
            as text "cluster bootstrap, `reps' reps, `level'% CI)"
        matlist `attm', border(rows) format(%9.4f)
        di as text _n "Average causal response ACRT(d)"
        matlist `acrm', border(rows) format(%9.4f)
        return matrix attd = `attm'
        return matrix acrt = `acrm'
        return scalar degree = `degree'
        return scalar dmin = `dmin'
        return scalar dmax = `dmax'
        restore
        exit
    }

    * ---------- first difference dy = y(post) - y(pre) ----------
    tempvar dy
    sort `unit' `time'
    by `unit': gen double `dy' = `depvar'[_N] - `depvar'[1]
    qui by `unit': keep if _n == _N        // one row per unit

    * ---------- dose grid ----------
    qui summarize `dose' if `dose' > 0, meanonly
    local dmin = r(min)
    local dmax = r(max)

    * ---------- explicit knots ----------
    local knotmatname ""
    if "`knots'" != "" {
        matrix _knots = (`knots')
        local knotmatname "_knots"
    }

    * ---------- point estimate (Mata) ----------
    tempname attm acrm
    matrix `attm' = J(`npoints', 5, .)
    matrix `acrm' = J(`npoints', 5, .)
    mata: _contdid_run("`dose'", "`dy'", `degree', "`knotmatname'", `nknots', `npoints', `dmin', `dmax', "_att", "_acrt")
    forvalues k = 1/`npoints' {
        local dk = `dmin' + (`dmax' - `dmin') * (`k' - 1) / (`npoints' - 1)
        matrix `attm'[`k',1] = `dk'
        matrix `attm'[`k',2] = _att[1,`k']
        matrix `acrm'[`k',1] = `dk'
        matrix `acrm'[`k',2] = _acrt[1,`k']
    }

    * ---------- uniform confidence band (cband) ----------
    if "`cband'" != "" {
        set seed `seed'
        mata: _contdid_cband("`dose'", "`dy'", `degree', "`knotmatname'", `nknots', `npoints', `dmin', `dmax', `reps', `level')
    }

    * ---------- cluster bootstrap ----------
    if `reps' > 0 {
        set seed `seed'
        if "`cluster'" != "" {
            local bscmd "bsample, cluster(`cluster')"
        }
        else {
            local bscmd "bsample"
        }
        tempfile est
        qui save `est'

        matrix boot = J(`reps', `npoints', .)
        matrix boota = J(`reps', `npoints', .)
        forvalues b = 1/`reps' {
            qui use `est', clear
            qui `bscmd'
            mata: _contdid_run("`dose'", "`dy'", `degree', "`knotmatname'", `nknots', `npoints', `dmin', `dmax', "_batt", "_bacrt")
            forvalues k = 1/`npoints' {
                matrix boot[`b',`k']  = _batt[1,`k']
                matrix boota[`b',`k'] = _bacrt[1,`k']
            }
        }

        qui use `est', clear
        local plo = (100 - `level') / 2
        local phi = 100 - `plo'
        svmat boot, names(b_)
        forvalues k = 1/`npoints' {
            qui summarize b_`k'
            matrix `attm'[`k',3] = r(sd)
            qui centile b_`k', centile(`plo' `phi')
            matrix `attm'[`k',4] = r(c_1)
            matrix `attm'[`k',5] = r(c_2)
        }
        svmat boota, names(c_)
        forvalues k = 1/`npoints' {
            qui summarize c_`k'
            matrix `acrm'[`k',3] = r(sd)
            qui centile c_`k', centile(`plo' `phi')
            matrix `acrm'[`k',4] = r(c_1)
            matrix `acrm'[`k',5] = r(c_2)
        }
    }

    matrix colnames `attm' = d ATT se lb ub
    matrix colnames `acrm' = d ACRT se lb ub

    * ---------- display ----------
    di as text _n "Dose-response ATT(d)  (B-spline degree `degree'; " ///
        as text "cluster bootstrap, `reps' reps, `level'% CI)"
    matlist `attm', border(rows) format(%9.4f)
    di as text _n "Average causal response ACRT(d)  (derivative of ATT(d))"
    matlist `acrm', border(rows) format(%9.4f)

    if "`cband'" != "" {
        di as text _n "Uniform confidence band (sup-t): crit_ATT = " ///
            as result %9.3f scalar(crit_att) as text ",  crit_ACRT = " ///
            as result %9.3f scalar(crit_acrt)
        di as text "ATT(d) with uniform band [cb_lb, cb_ub]:"
        matlist cb_att, border(rows) format(%9.4f)
    }

    * ---------- graph ----------
    if "`graph'" != "" {
        qui clear
        qui set obs `npoints'
        gen double d   = .
        gen double att = .
        gen double lb  = .
        gen double ub  = .
        gen double acrt = .
        forvalues k = 1/`npoints' {
            qui replace d    = `attm'[`k',1] in `k'
            qui replace att  = `attm'[`k',2] in `k'
            qui replace lb   = `attm'[`k',4] in `k'
            qui replace ub   = `attm'[`k',5] in `k'
            qui replace acrt = `acrm'[`k',2] in `k'
        }
        twoway (rarea ub lb d, color(gs13)) ///
               (line att d, lcolor(navy) lwidth(medthick)) ///
               (line acrt d, lcolor(maroon) lwidth(medthick) yaxis(2)), ///
            legend(order(2 "ATT(d)" 3 "ACRT(d)") rows(1)) ///
            title("Dose-response: ATT(d) and ACRT(d)") ///
            xtitle("Dose (d)") ytitle("ATT(d)", axis(1)) ytitle("ACRT(d)", axis(2)) ///
            note("B-spline degree `degree'; cluster bootstrap `level'% CI", size(small))
    }

    * ---------- returns ----------
    return matrix attd = `attm'
    return matrix acrt = `acrm'
    return scalar degree = `degree'
    return scalar dmin   = `dmin'
    return scalar dmax   = `dmax'
    if "`cband'" != "" {
        return scalar crit_att = crit_att
        return scalar crit_acrt = crit_acrt
        return matrix cb_att  = cb_att
        return matrix cb_acrt = cb_acrt
    }

    restore
end
