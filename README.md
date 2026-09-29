# contdid

**Difference-in-Differences with a Continuous Treatment for Stata**

[![Stata 16+](https://img.shields.io/badge/Stata-16%2B-blue.svg)](https://www.stata.com/)
[![License: AGPL-3.0](https://img.shields.io/badge/License-AGPL--3.0-blue.svg)](LICENSE)
[![Version: 0.1.0](https://img.shields.io/badge/Version-0.1.0-green.svg)](CHANGELOG.md)

`contdid` implements the **difference-in-differences (DiD) estimator with a continuous
treatment** of Callaway, Goodman-Bacon & Sant'Anna (2024/2025), filling a gap in the Stata
ecosystem: the authors' reference implementation (`contdid`) is available only in **R**.

## Overview

In many DiD designs the treatment does not simply "turn on" — it has a *dose* or *intensity*
(policy generosity, subsidy size, pollution exposure, financial-development index, …).
`contdid` estimates the **dose-response function** ATT(d): the average treatment effect on the
treated as a function of the continuous treatment dose.

Under parallel trends, the dose-specific effect is identified as

> **ATT(d) = E[ΔY | D = d] − E[ΔY | D = 0]**,

where ΔY is the pre-to-post change in the outcome, D is the continuous dose, and D = 0 marks
untreated units. The current version (v1) estimates this with a **linear-in-dose** specification
(the default of the reference `contdid` package, `num_knots = 0, degree = 1`) and reports
pointwise confidence intervals via **cluster bootstrap**.

> **v1 scope:** two-period panel / repeated cross-section, continuous dose, ATT(d) only.
> Staggered adoption, ACRT (causal-response/slope) parameters, uniform confidence bands and
> the CCK data-driven sieve estimator are on the roadmap (see [CHANGELOG](CHANGELOG.md)).

## Installation

```stata
net install contdid, from("https://raw.githubusercontent.com/niuniuhaoyu/contdid/main/") replace
```

Requires Stata 16 or later. No additional dependencies.

## Quick start

```stata
* Load a two-period panel with variables id, t (0/1), y (outcome), d (continuous dose ≥ 0)
use contdid_sim.dta, clear

* Estimate the dose-response ATT(d) with a plot
contdid y, unit(id) time(t) dose(d) npoints(20) reps(999) seed(12345) graph
```

| Option | Default | Description |
|---|---|---|
| `unit(varname)` | — | unit / panel identifier (required) |
| `time(varname)` | — | time variable, two periods (required) |
| `dose(varname)` | — | continuous treatment dose, ≥ 0, 0 = untreated (required) |
| `npoints(#)` | 20 | number of evaluation doses |
| `reps(#)` | 999 | cluster-bootstrap replications |
| `seed(#)` | 12345 | random seed |
| `cluster(varname)` | `unit()` | cluster variable for bootstrap |
| `level(#)` | 95 | confidence level (%) |
| `graph` | off | plot the dose-response curve with CI band |

## Method

- **Estimator** (linear-in-dose): fit ΔY = a + b·D by OLS on treated units (D > 0), then
  ATT(d) = (a + b·d) − mean(ΔY | D = 0). This matches the reference `contdid` default.
- **Inference**: cluster bootstrap (default cluster = unit), percentile pointwise confidence
  intervals.
- **Interpretation caveat**: ATT(d) is the *local* effect for dose group d. The *global*
  dose-response ATT(d) requires the stronger "strong parallel trends" assumption; average
  causal responses (ACRT) are not yet implemented. See the paper for details.

## Reproducibility

- `examples/contdid_simdata.do` — generates the simulated dataset (`data/contdid_sim.dta`)
  whose true dose-response is ATT(d) = 0.5·d.
- `examples/contdid_example.do` — one-click example.
- `examples/reference/run_contdid_check.R` — independent base-R reimplementation of the same
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
  version = {0.1.0},
  url     = {https://github.com/niuniuhaoyu/contdid}
}
```

## License

[AGPL-3.0](LICENSE)
