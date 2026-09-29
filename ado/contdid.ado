*! contdid: Difference-in-Differences with a Continuous Treatment
*! version 0.1.0  2026-09-29  Haoyu Niu
*! Implements Callaway, Goodman-Bacon & Sant'Anna (2024/2025)
*! Two-period continuous-treatment DiD, linear-in-dose ATT(d)
*!
*! ATT(d) = E[ΔY | D=d] - E[ΔY | D=0]
*!   linear fit: ΔY = a + b*D  on treated (D>0) units
*!   baseline:   mean(ΔY | D=0)  (untreated units)

program define contdid, rclass
    version 16

    syntax varlist(max=1 numeric) [if] [in], ///
        unit(varname numeric) ///           individual id
        time(varname numeric) ///           time period (two periods)
        dose(varname numeric) ///           continuous treatment dose
        [npoints(integer 20) ///             number of evaluation doses
         level(real 95) ///                  confidence level (%)
         seed(integer 12345) ///             RNG seed (used in bootstrap)
         reps(integer 999) ///               bootstrap replications
         cluster(varname) ///                cluster variable (default = unit)
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
    if `ntvals' != 2 {
        di as error "contdid: v1 supports only two time periods; found `ntvals' distinct values"
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

    * ---------- first difference dy = y(post) - y(pre) ----------
    tempvar dy
    sort `unit' `time'
    by `unit': gen double `dy' = `depvar'[_N] - `depvar'[1]
    qui by `unit': keep if _n == _N        // one row per unit

    * ---------- linear-in-dose among treated ----------
    qui regress `dy' `dose' if `dose' > 0
    local b_d    = _b[`dose']
    local b_cons = _b[_cons]

    * ---------- untreated baseline ----------
    qui summarize `dy' if `dose' == 0, meanonly
    local m0 = r(mean)

    * ---------- dose grid and ATT(d) ----------
    qui summarize `dose' if `dose' > 0, meanonly
    local dmin = r(min)
    local dmax = r(max)

    tempname attm
    matrix `attm' = J(`npoints', 5, .)
    forvalues k = 1/`npoints' {
        local dk  = `dmin' + (`dmax' - `dmin') * (`k' - 1) / (`npoints' - 1)
        local att = `b_cons' + `b_d' * `dk' - `m0'
        matrix `attm'[`k',1] = `dk'
        matrix `attm'[`k',2] = `att'
        // cols 3-5 (se, lb, ub) filled by cluster bootstrap (Task 4)
    }
    matrix colnames `attm' = d ATT se lb ub

    di as text _n "Dose-response ATT(d)  (linear-in-dose; SE/CI in Task 4)"
    matlist `attm', border(rows) format(%9.4f)

    * ---------- returns ----------
    return matrix attd = `attm'
    return scalar b_d    = `b_d'
    return scalar b_cons = `b_cons'
    return scalar m0     = `m0'

    restore
end
