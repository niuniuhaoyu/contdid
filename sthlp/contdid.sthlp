{smcl}
{* 29 Sep 2026}{...}
{hline}
{p 4 8 2}{bf:contdid} — Difference-in-Differences with a continuous treatment{right:version 0.2.0}
{hline}

{title:Title}

{p 4 4 2}
{cmd:contdid} — Estimate the dose-response function ATT(d) and the average
causal response ACRT(d) for a difference-in-differences design with a
continuous treatment.

{title:Syntax}

{p 8 12 2}
{cmd:contdid} {it:depvar} {ifin}, {cmdab:unit:(}{it:varname}{cmd:)}
{cmdab:time:(}{it:varname}{cmd:)} {cmdab:dose:(}{it:varname}{cmd:)}
{cmd:[}{cmdab:gvar:(}{it:varname}{cmd:)} {cmdab:degree:(}{it:#}{cmd:)}
{cmdab:nknots:(}{it:#}{cmd:)}
{cmdab:knots:(}{it:numlist}{cmd:)} {cmdab:npoints:(}{it:#}{cmd:)}
{cmdab:level:(}{it:#}{cmd:)} {cmdab:reps:(}{it:#}{cmd:)}
{cmdab:seed:(}{it:#}{cmd:)} {cmdab:cluster:(}{it:varname}{cmd:)}
{cmd:cband} {cmd:graph}{cmd:]}

{title:Description}

{p 4 4 2}
{cmd:contdid} implements the two-period difference-in-differences estimator
with a continuous treatment of Callaway, Goodman-Bacon, and Sant'Anna
(2024, NBER WP 32117). Under parallel trends the dose-specific average
treatment effect on the treated is identified as

{p 12 12 2}
ATT(d) = E[ΔY | D = d] − E[ΔY | D = 0]

{p 4 4 2}
where ΔY is the pre-to-post change in {it:depvar}, D is the continuous dose
({cmd:dose()}), and D = 0 marks untreated units. {cmd:contdid} models the
dose-response with a B-spline (OLS of ΔY on a B-spline basis of D among
treated units, minus the untreated-group mean of ΔY), estimates ATT(d) and its
derivative ACRT(d), and reports pointwise confidence intervals from a cluster
bootstrap; {cmd:cband} adds a uniform confidence band via a multiplier
bootstrap.

{title:Options}

{p 4 4 2}{cmdab:unit:(}{it:varname}{cmd:)} identifies the panel units.
{p 4 4 2}{cmdab:time:(}{it:varname}{cmd:)} is the time variable; it must take
exactly two distinct values.
{p 4 4 2}{cmdab:dose:(}{it:varname}{cmd:)} is the continuous treatment dose.
It must be non-negative, with 0 indicating untreated units; a group with
{cmd:dose} = 0 is required.
{p 4 4 2}{cmdab:gvar:(}{it:varname}{cmd:)} is the treatment timing (the period
when a unit is first treated; 0 = never treated). Specifying {cmd:gvar()} enables
staggered adoption with multiple treatment cohorts.
{p 4 4 2}{cmdab:degree:(}{it:#}{cmd:)} sets the B-spline degree (default 1 =
linear; 2 = quadratic; 3 = cubic).
{p 4 4 2}{cmdab:nknots:(}{it:#}{cmd:)} sets the number of interior knots,
placed at dose quantiles (default 0 = global polynomial).
{p 4 4 2}{cmdab:knots:(}{it:numlist}{cmd:)} specifies interior knots
explicitly (overrides {cmd:nknots()}).
{p 4 4 2}{cmdab:npoints:(}{it:#}{cmd:)} sets the number of evaluation doses
(default 20, minimum 2).
{p 4 4 2}{cmdab:level:(}{it:#}{cmd:)} sets the confidence level in percent
(default 95).
{p 4 4 2}{cmdab:reps:(}{it:#}{cmd:)} sets the number of bootstrap
replications (default 999); also used for the {cmd:cband} multiplier draws.
{p 4 4 2}{cmdab:seed:(}{it:#}{cmd:)} sets the random seed (default 12345).
{p 4 4 2}{cmdab:cluster:(}{it:varname}{cmd:)} resamples clusters instead of
units (default is the unit identifier).
{p 4 4 2}{cmd:cband} reports a uniform confidence band (sup-t) via a
multiplier bootstrap.
{p 4 4 2}{cmd:graph} plots ATT(d) and ACRT(d).

{title:Saved results}

{p 4 4 2}
{cmd:r(attd)} and {cmd:r(acrt)} are {it:npoints} × 5 matrices with columns
{cmd:d}, {cmd:ATT}/{cmd:ACRT}, {cmd:se}, {cmd:lb}, {cmd:ub}.
{cmd:r(degree)}, {cmd:r(dmin)}, and {cmd:r(dmax)} return the spline degree and
the dose range.
{p 4 4 2}
With {cmd:cband}: {cmd:r(crit_att)} and {cmd:r(crit_acrt)} are the sup-t
critical values; {cmd:r(cb_att)} and {cmd:r(cb_acrt)} are {it:npoints} × 4
matrices with columns {cmd:d}, estimate, {cmd:cb_lb}, {cmd:cb_ub}.

{title:Examples}

{p 4 4 2}
Linear dose-response of {cmd:y} on continuous dose {cmd:d}, with a plot:

{phang2}{cmd:. contdid y, unit(id) time(t) dose(d) npoints(20) reps(999) seed(1) graph}{p_end}

{p 4 4 2}
Cubic B-spline with 2 interior knots and a uniform confidence band:

{phang2}{cmd:. contdid y, unit(id) time(t) dose(d) degree(3) nknots(2) reps(999) seed(1) cband}{p_end}

{title:References}

{p 4 4 2}
Callaway, B., A. Goodman-Bacon, and P. H. C. Sant'Anna. 2024.
Difference-in-Differences with a Continuous Treatment. NBER Working Paper 32117.

{p 4 4 2}
Niu, H. 2026. contdid: Difference-in-Differences with a Continuous Treatment
for Stata. Version 0.2.0.

{title:Also see}

{p 4 4 2}
Help: {help csdid} (binary staggered DiD), {help drdid} (doubly robust DiD).
