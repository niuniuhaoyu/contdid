# Changelog

All notable changes to this project will be documented in this file.

## [1.0.0] - 2026-09-29

First stable release. Two-period and staggered continuous-treatment DiD, B-spline
dose-response ATT(d) and ACRT(d), cluster-bootstrap pointwise intervals, and a
multiplier-bootstrap uniform confidence band.

## [0.3.0] - 2026-09-29

### Added
- Staggered adoption support (`gvar()` option): group-time dose-response estimation across
  multiple treatment cohorts, aggregated into ATT(d) and ACRT(d), with cluster bootstrap.

## [0.2.0] - 2026-09-29

### Added
- B-spline dose-response modeling (`degree()`, `knots()`, `nknots()` options); default `degree(1)` preserves v0.1.0 behavior.
- ACRT(d) — average causal response (derivative of ATT(d)), returned as `r(acrt)`.
- Uniform confidence band (`cband` option) via influence-function multiplier bootstrap (sup-t).
- Mata B-spline basis/derivative (Cox-de Boor), verified against R `splines2`.

## [0.1.0] - 2026-09-29

### Added
- `contdid` command: difference-in-differences with a continuous treatment (two-period, linear-in-dose).
- Dose-response estimation of ATT(d) via linear-in-dose fit with untreated-group baseline.
- Cluster bootstrap pointwise confidence intervals.
- Dose-response plot (`graph` option).
- Simulated example dataset (`contdid_sim.dta`).
