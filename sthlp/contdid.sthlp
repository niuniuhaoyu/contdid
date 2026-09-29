{smcl}
{* 29 Sep 2026}{...}
{hline}
{p 4 8 2}{bf:contdid} — Difference-in-Differences with a continuous treatment{right:version 0.1.0}
{hline}

{title:Title}

{p 4 4 2}
{cmd:contdid} — Estimate the dose-response function ATT(d) for a
difference-in-differences design with a continuous treatment.

{title:Syntax}

{p 8 12 2}
{cmd:contdid} {it:depvar} {ifin}, {cmdab:unit:(}{it:varname}{cmd:)}
{cmdab:time:(}{it:varname}{cmd:)} {cmdab:dose:(}{it:varname}{cmd:)}
{cmd:[}{cmdab:npoints:(}{it:#}{cmd:)} {cmdab:level:(}{it:#}{cmd:)}
{cmdab:reps:(}{it:#}{cmd:)} {cmdab:seed:(}{it:#}{cmd:)}
{cmdab:cluster:(}{it:varname}{cmd:)} {cmd:graph}{cmd:]}

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
({cmd:dose()}), and D = 0 marks untreated units. {cmd:contdid} estimates this
with a linear-in-dose specification (OLS of ΔY on D among treated units,
minus the untreated-group mean of ΔY) and reports pointwise confidence
intervals from a cluster bootstrap.

{title:Options}

{p 4 4 2}{cmdab:unit:(}{it:varname}{cmd:)} identifies the panel units.
{p 4 4 2}{cmdab:time:(}{it:varname}{cmd:)} is the time variable; it must take
exactly two distinct values (v1 supports two periods only).
{p 4 4 2}{cmdab:dose:(}{it:varname}{cmd:)} is the continuous treatment dose.
It must be non-negative, with 0 indicating untreated units; a group with
{cmd:dose} = 0 is required.
{p 4 4 2}{cmdab:npoints:(}{it:#}{cmd:)} sets the number of evaluation doses
(default 20, minimum 2).
{p 4 4 2}{cmdab:level:(}{it:#}{cmd:)} sets the confidence level in percent
(default 95).
{p 4 4 2}{cmdab:reps:(}{it:#}{cmd:)} sets the number of bootstrap
replications (default 999).
{p 4 4 2}{cmdab:seed:(}{it:#}{cmd:)} sets the random seed (default 12345).
{p 4 4 2}{cmdab:cluster:(}{it:varname}{cmd:)} resamples clusters instead of
units (default is the unit identifier).
{p 4 4 2}{cmd:graph} plots the dose-response curve with its confidence band.

{title:Saved results}

{p 4 4 2}
{cmd:r(attd)} is an {it:npoints} × 5 matrix with columns {cmd:d}, {cmd:ATT},
{cmd:se}, {cmd:lb}, {cmd:ub}.
{cmd:r(b_d)}, {cmd:r(b_cons)}, and {cmd:r(m0)} return the dose coefficient,
the intercept, and the untreated-group baseline, respectively.

{title:Examples}

{p 4 4 2}
Estimate the dose-response of {cmd:y} on continuous dose {cmd:d}:

{phang2}{cmd:. contdid y, unit(id) time(t) dose(d) npoints(20) reps(999) seed(1) graph}{p_end}

{title:References}

{p 4 4 2}
Callaway, B., A. Goodman-Bacon, and P. H. C. Sant'Anna. 2024.
Difference-in-Differences with a Continuous Treatment. NBER Working Paper 32117.

{p 4 4 2}
Niu, H. 2026. contdid: Difference-in-Differences with a Continuous Treatment
for Stata. Version 0.1.0.

{title:Also see}

{p 4 4 2}
Help: {help csdid} (binary staggered DiD), {help drdid} (doubly robust DiD).
