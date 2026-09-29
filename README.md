# contdid

**Difference-in-Differences with a Continuous Treatment for Stata**

[![Stata 16+](https://img.shields.io/badge/Stata-16%2B-blue.svg)](https://www.stata.com/)
[![License: AGPL-3.0](https://img.shields.io/badge/License-AGPL--3.0-blue.svg)](LICENSE)
[![Version: 0.3.0](https://img.shields.io/badge/Version-0.3.0-green.svg)](CHANGELOG.md)

`contdid` implements the **difference-in-differences (DiD) estimator with a continuous
treatment** of Callaway, Goodman-Bacon & Sant'Anna (2024/2025), filling a gap in the Stata
ecosystem: the authors' reference implementation (`contdid`) is available only in **R**.

## Overview

In many DiD designs the treatment does not simply "turn on" — it has a *dose* or *intensity*
(policy generosity, subsidy size, pollution exposure, financial-development index, …).
`contdid` estimates the **dose-response function** ATT(d) and its derivative, the **average
causal response** ACRT(d).

Under parallel trends, the dose-specific effect is identified as

> **ATT(d) = E[ΔY | D = d] − E[ΔY | D = 0]**,

where ΔY is the pre-to-post change in the outcome, D is the continuous dose, and D = 0 marks
untreated units. The dose-response is modeled with a **B-spline** (default: linear), and
ACRT(d) is the derivative of the fitted curve. Pointwise confidence intervals come from a
cluster bootstrap, and an optional **uniform confidence band** (`cband`) covers the whole
curve simultaneously via a multiplier bootstrap.

> **Scope:** two-period panel / repeated cross-section, or **staggered adoption** (with `gvar()`),
> continuous dose, ATT(d) and ACRT(d). Covariates and the CCK data-driven sieve estimator are on
> the roadmap (see [CHANGELOG](CHANGELOG.md)).

## Installation

```stata
net install contdid, from("https://raw.githubusercontent.com/niuniuhaoyu/contdid/main/") replace
```

Requires Stata 16 or later. No additional dependencies.

## Quick start

```stata
* Load a two-period panel with variables id, t (0/1), y (outcome), d (continuous dose ≥ 0)
use contdid_sim.dta, clear

* Linear dose-response (default), with plot
contdid y, unit(id) time(t) dose(d) npoints(20) reps(999) seed(12345) graph

* Cubic B-spline with 2 interior knots, plus a uniform confidence band
contdid y, unit(id) time(t) dose(d) degree(3) nknots(2) reps(999) seed(12345) cband graph

* Staggered adoption: units treated at different times (g = treatment period, 0 = never)
contdid y, unit(id) time(t) dose(d) gvar(g) npoints(20) reps(999) seed(12345)
```

| Option | Default | Description |
|---|---|---|
| `unit(varname)` | — | unit / panel identifier (required) |
| `time(varname)` | — | time variable, two periods (required) |
| `dose(varname)` | — | continuous treatment dose, ≥ 0, 0 = untreated (required) |
| `gvar(varname)` | — | treatment timing (0 = never treated); enables staggered adoption |
| `degree(#)` | 1 | B-spline degree (1 = linear, 2 = quadratic, 3 = cubic) |
| `nknots(#)` | 0 | number of interior knots (placed at dose quantiles) |
| `knots(numlist)` | — | explicit interior knots (overrides `nknots`) |
| `npoints(#)` | 20 | number of evaluation doses |
| `reps(#)` | 999 | cluster-bootstrap replications (and multiplier draws) |
| `seed(#)` | 12345 | random seed |
| `cluster(varname)` | `unit()` | cluster variable for bootstrap |
| `level(#)` | 95 | confidence level (%) |
| `cband` | off | uniform confidence band (sup-t, multiplier bootstrap) |
| `graph` | off | plot ATT(d) and ACRT(d) |

## Method

- **Estimator**: fit ΔY on a B-spline basis of the dose (OLS, treated units), then
  ATT(d) = fitted(d) − mean(ΔY | D = 0). ACRT(d) is the derivative of the fitted curve.
  With `degree(1)` and `nknots(0)` this reduces to the linear-in-dose estimator of v0.1.0.
- **Pointwise inference**: cluster bootstrap (default cluster = unit), percentile confidence
  intervals.
- **Uniform inference** (`cband`): influence-function-based multiplier bootstrap (sup-t) gives
  a band that covers the whole dose-response curve simultaneously.
- **Interpretation caveat**: ATT(d) is the *local* effect for dose group d; the *global*
  dose-response requires the stronger "strong parallel trends" assumption. ACRT(d) is a
  derivative and is estimated with less precision than ATT(d). See the paper for details.

## Reproducibility

- `examples/contdid_simdata.do` — generates the simulated dataset (`data/contdid_sim.dta`)
  whose true dose-response is ATT(d) = 0.5·d.
- `examples/contdid_example.do` — one-click example.
- `examples/reference/dump_splines2.R` — dumps R `splines2` B-spline basis/derivative as the
  golden reference for the Mata implementation.
- `examples/reference/run_contdid_check.R` — independent base-R reimplementation of the linear
  estimator; matches the Stata output to 6 decimals.

## Citation

Method:

```bibtex
@article{callaway2024continuous,
  title   = {Difference-in-Differences with a Continuous Treatment},
  author  = {Callaway, Brantly and Goodman-Bacon, Andrew and Sant'Anna, Pedro H. C.},
  year    = {2024},
  note    = {NBER Working Paper 32117},
  url     = {https://doi.org/10.3386/w32117}
}
```

Software:

```bibtex
@software{niu2026contdid,
  title   = {contdid: Difference-in-Differences with a Continuous Treatment for Stata},
  author  = {Haoyu Niu},
  year    = {2026},
  version = {0.2.0},
  url     = {https://github.com/niuniuhaoyu/contdid}
}
```

## License

[AGPL-3.0](LICENSE)
