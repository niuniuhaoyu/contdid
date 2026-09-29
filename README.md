# contdid

**Difference-in-Differences with a Continuous Treatment for Stata**

[![Stata 16+](https://img.shields.io/badge/Stata-16%2B-blue.svg)](https://www.stata.com/)
[![License: AGPL-3.0](https://img.shields.io/badge/License-AGPL--3.0-blue.svg)](LICENSE)
[![Version: 1.0.0](https://img.shields.io/badge/Version-1.0.0-green.svg)](CHANGELOG.md)

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

## Methods

### Setup and notation

Consider a panel of units *i* = 1, …, *n* observed over periods *t* = 1, …, *T*. Each unit has a
**dose** *Dᵢ* ≥ 0 (constant over time; *D* = 0 marks untreated units) and, in the staggered
design, a **treatment timing** *Gᵢ* (the first period in which the unit is treated; *G* = 0 if
never treated). Let *Yᵢₜ*(*d*) be the potential outcome under dose *d*. Under no anticipation
the observed outcomes satisfy *Yᵢₜ* = *Yᵢₜ*(0) before treatment and *Yᵢₜ* = *Yᵢₜ*(*Dᵢ*)
afterwards. Write Δ*Y* = *Y*₂ − *Y*₁ for the pre-to-post change.

### Identification

The key identifying assumption is **parallel trends**:

> For all *d* > 0,  E[*Y*₂(0) − *Y*₁(0) | *D* = *d*] = E[*Y*₂(0) − *Y*₁(0) | *D* = 0].

Under this assumption the **dose-specific average treatment effect on the treated** is
identified (Theorem 3.1 of the paper) as

> **ATT(*d* | *d*) = E[Δ*Y* | *D* = *d*] − E[Δ*Y* | *D* = 0]**.

Here ATT(*d* | *d*) is the *local* effect for dose group *d*. The *global* dose-response
ATT(*d*) = E[*Y*₂(*d*) − *Y*₂(0) | *D* > 0] requires the stronger **strong parallel trends**
assumption, which rules out selection into dose levels. The **average causal response**
ACRT(*d*) is the derivative of the dose-response with respect to *d*.

### Estimation

`contdid` models the dose-response nonparametrically with a **B-spline**. Among treated units
(*D* > 0) it fits

> Δ*Y* = *α* + Σₖ βₖ *B*ₖ(*D*) + *ε*,

where *B*ₖ are B-spline basis functions (degree and interior knots chosen by the user). The
estimators are

> ATT(*d*) = *α̂* + Σₖ β̂ₖ *B*ₖ(*d*) − *m̂*₀,   with *m̂*₀ = mean(Δ*Y* | *D* = 0),
>
> ACRT(*d*) = Σₖ β̂ₖ *B*′ₖ(*d*).

With `degree(1)` and `nknots(0)` the basis is linear, and this reduces to the simple
linear-in-dose estimator E[Δ*Y* | *D* = *d*] − E[Δ*Y* | *D* = 0].

### Inference

- **Pointwise confidence intervals** come from a cluster bootstrap (default cluster = unit):
  resample units with replacement, re-estimate, and report the percentile interval of the
  bootstrap estimates.
- **Uniform confidence band** (`cband`) covers the whole dose-response curve simultaneously.
  It is built from the influence function of the spline estimator with a multiplier bootstrap
  (Rademacher weights); the sup-t critical value is the (1 − *α*) quantile of the studentized
  max over the evaluation grid.

### Staggered adoption

With `gvar()`, units are first treated at potentially different times. For each treatment
cohort *g* and post-treatment period *t*, `contdid` forms a "2×2"-style comparison of cohort
*g* against the not-yet-treated (units with *G* > *t* or *G* = 0), differences outcomes
Δ*Y* = *Yₜ* − *Y*_{*g*−1}, and estimates the cohort-time dose-response ATT(*g*, *t*; *d*) with
the same B-spline estimator. The reported ATT(*d*) and ACRT(*d*) aggregate these cohort-time
effects across all (*g*, *t*) pairs, weighted by cohort size.

### Interpretation caveats

- ATT(*d*) is a *local* dose effect; a *global* interpretation requires strong parallel trends.
- ACRT(*d*) is a derivative and is estimated with less precision than ATT(*d*).
- The dose *D* must be non-negative and constant within unit; a group with *D* = 0 is required.
- As with any DiD design, credibility rests on the plausibility of parallel trends
  (flat pre-treatment trends).

## Usage guide

### A minimal two-period example

```stata
* Two-period panel: id, t (0/1), y, d (continuous dose)
use contdid_sim.dta, clear
contdid y, unit(id) time(t) dose(d)
```

`contdid` reports two tables: **ATT(d)** (the dose-response) and **ACRT(d)** (its derivative).
The dose-response is the effect of receiving dose *d* versus dose 0; ACRT(d) is the effect of a
one-unit increase in the dose at level *d*. With `graph`, both are plotted; with `cband`, the
band covers the curve simultaneously rather than point-by-point.

### Choosing the spline

- `degree(1)` (default) imposes a linear dose-response — the most parsimonious, and a good
  starting point.
- `degree(2)`/`degree(3)` with `nknots(1)`–`nknots(3)` lets the data reveal curvature; compare
  the fit visually and report the specification you rely on.
- `knots(0.3 0.6)` fixes the knots at specific dose values for reproducibility across
  specifications.

### Staggered adoption

```stata
* Staggered: g = treatment period (0 = never treated)
contdid y, unit(id) time(t) dose(d) gvar(g) cband graph
```

### Reproducing the published figure

```stata
do examples/make_figure.do   // generates examples/contdid_figure.png
```

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
  version = {1.0.0},
  url     = {https://github.com/niuniuhaoyu/contdid}
}
```

## License

[AGPL-3.0](LICENSE)
