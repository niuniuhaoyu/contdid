# Changelog

All notable changes to this project will be documented in this file.

## [1.1.0] - 2026-10-07

### Added (v3a)
- `covariates(varlist)` option: conditional (strong) parallel trends estimation
  of ATT(d)/ACRT(d), following Callaway-Goodman-Bacon-Sant'Anna (2024/2025),
  Supplemental Appendix SI.3 (Assumption SPT-X, Proposition S3). Joint sieve
  with dose-by-covariate interactions; `ATT(d) = E_X[ATT_x(d) | D>0]`.
  Supports `cband` (influence-function multiplier bootstrap that accounts for
  both the OLS parameters and the treated covariate means).
  Verified on a conditional-parallel-trends DGP (`examples/_test_v3a_cov.do`)
  against an independent R prototype (`examples/reference/v3a_prototype.R`).
- Not yet supported (planned): `covariates()` with staggered adoption (`gvar()`).

### Added (v3b, scoped)
- `dose_est_method(dds)`: data-driven sieve — selects the number of interior knots
  (0..`maxknots()`) by leave-one-out CV on the treated fit, then estimates ATT(d)/ACRT(d).
  NOTE: this is a **tractable data-driven sieve, not the exact Chen-Christensen-Kankanala /
  npiv "cck" estimator** (which uses Tikhonov regularization and a data-driven J). Exact
  parity with R `contdid`'s `dose_est_method="cck"` is **not** achieved. Verified on a
  nonlinear DGP (`examples/_test_v3b_dds.do`).

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
