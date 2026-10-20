{smcl}
{* *! version 2.0.0  17oct2026}{...}
{* *! Mata port of the Python reference implementation (ivmediate 2.0.0)}{...}
{viewerjumpto "Syntax" "ivmediate##syntax"}{...}
{viewerjumpto "Description" "ivmediate##description"}{...}
{viewerjumpto "Where to find the theory" "ivmediate##papermap"}{...}
{viewerjumpto "What's new in 2.0" "ivmediate##whatsnew"}{...}
{viewerjumpto "Options" "ivmediate##options"}{...}
{viewerjumpto "Remarks" "ivmediate##remarks"}{...}
{viewerjumpto "Conventions" "ivmediate##conventions"}{...}
{viewerjumpto "Examples" "ivmediate##examples"}{...}
{viewerjumpto "Stored results" "ivmediate##results"}{...}
{viewerjumpto "Citation" "ivmediate##citation"}{...}
{viewerjumpto "References" "ivmediate##references"}{...}
{viewerjumpto "Authors" "ivmediate##authors"}{...}
{viewerjumpto "License" "ivmediate##license"}{...}
{cmd:help ivmediate}{right: (version 2.0.0)}
{hline}

{title:Title}

{p2colset 5 18 20 2}{...}
{p2col :{cmd:ivmediate} {hline 2}}Causal mediation analysis for linear
instrumental-variables models (no separate instrument for the mediator is
required){p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:ivmediate} {depvar} [{indepvars}] {ifin} {weight}{cmd:,}
{cmdab:med:iator}{cmd:(}{it:{help varname}}{cmd:)}
{cmdab:treat:ment}{cmd:(}{it:{help varname}}{cmd:)}
{cmdab:inst:rument}{cmd:(}{it:{help varlist}}{cmd:)}
[{it:options}]

{synoptset 30 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{p2coldent:* {opt med:iator(varname)}}single mediator variable{p_end}
{p2coldent:* {opt treat:ment(varname)}}single treatment variable{p_end}
{p2coldent:* {opt inst:rument(varlist)}}one {it:or more} excluded instruments{p_end}
{synopt: {opt a:bsorb(varlist)}}one or two categorical variables (numeric or string) to absorb as fixed effects{p_end}

{syntab:Standard errors}
{synopt: {cmd:vce(}{it:vcetype}{cmd:)}}{cmd:robust} (the default), {cmd:unadjusted}, {cmd:cluster} {it:clustvar}, or {cmd:conley}{p_end}
{synopt: {opt lat:itude(varname)}}latitude/y-coordinate; {cmd:vce(conley)} only{p_end}
{synopt: {opt lon:gitude(varname)}}longitude/x-coordinate; {cmd:vce(conley)} only{p_end}
{synopt: {opt cut:off(#)}}spatial distance cutoff; {cmd:vce(conley)} only{p_end}
{synopt: {opt kern:el(string)}}spatial kernel: {cmd:uniform} (default) or {cmd:bartlett}; {cmd:vce(conley)} only{p_end}
{synopt: {opt metr:ic(string)}}distance metric: {cmd:haversine} (default) or {cmd:euclidean}; {cmd:vce(conley)} only{p_end}

{syntab:Weak-instrument diagnostics and robust sets}
{synopt: {opt opt:au(#)}}bias fraction {it:tau} for the Olea-Pflueger critical value; default {cmd:0.10}{p_end}
{synopt: {opt opa:lpha(#)}}size {it:alpha} for the Olea-Pflueger critical value; default {cmd:0.05}{p_end}
{synopt: {opt noar}}skip the Anderson-Rubin confidence sets (computed by default){p_end}
{synopt: {opt ar}}accepted for compatibility; the sets are computed by default{p_end}
{synopt: {opt arp:oints(#)}}AR search-grid points (K>1); default {cmd:4001}{p_end}
{synopt: {opt arw:idth(#)}}AR search range in SE units; default {cmd:200}{p_end}

{syntab:Bootstrap}
{synopt: {opt b:ootstrap}}add bootstrap standard errors and confidence intervals{p_end}
{synopt: {opt bre:ps(#)}}bootstrap replications; default {cmd:999}{p_end}
{synopt: {opt bm:ethod(string)}}{cmd:block} (default with {cmd:vce(cluster)}) or {cmd:pairs} (default otherwise){p_end}
{synopt: {opt bci:type(string)}}{cmd:basic} (default), {cmd:percentile}, or {cmd:studentized}{p_end}
{synopt: {opt bse:ed(#)}}random-number seed; default {cmd:20260317}{p_end}
{synopt: {opt bootidx(matname)}}user-supplied resamples (reproducibility, cross-package checks){p_end}
{synopt: {opt bdiag:nostics}}also bootstrap RV(theta) and XRV(theta){p_end}

{syntab:Post-double-selection LASSO}
{synopt: {opt lasso}}select controls from a basis expansion by post-double-selection LASSO{p_end}
{synopt: {opt lpenalty(string)}}penalty rule; {cmd:plugin} (default and only supported value){p_end}
{synopt: {opt expandx(string)}}operators on the continuous controls; default {cmd:square cube}{p_end}
{synopt: {opt expanddum(string)}}operators on the factor dummies; default {cmd:interact_with_X}{p_end}
{synopt: {opt lkeep(patterns)}}expansion terms always kept (matched by name){p_end}
{synopt: {opt showsel}}display the per-equation selections (S_Y, S_T, S_M, S_Z){p_end}
{synopt: {opt lcslack(#)}}penalty slack {it:c} > 0; default {cmd:1.1}{p_end}
{synopt: {opt lalphap(#)}}penalty tail probability in (0, 1); default 0.1/ln({it:n}){p_end}
{synopt: {opt lniter(#)}}maximum penalty-loading iterations; default {cmd:15}{p_end}
{synopt: {opt ltol(#)}}convergence tolerance of the loading iterations; default {cmd:1e-6}{p_end}
{synopt: {opt llambda(#)}}fix the penalty level (skips the plug-in iterations){p_end}

{syntab:Reporting}
{synopt: {opt l:evel(#)}}confidence level; default {cmd:level(95)}{p_end}
{synopt: {opt nostar:s}}suppress significance stars in the estimates table{p_end}
{synopt: {opt starcut:offs(# # #)}}three p-value thresholds for stars; default {cmd:0.01 0.05 0.10}{p_end}
{synopt: {opt kappa:table}}print the kappa-adjusted estimates table{p_end}
{synopt: {opt r2:sensitivity}}print the R-squared-based sensitivity section{p_end}
{synopt: {opt resid:uals(stub)}}save the residuals as {it:stub}{cmd:_eT}, {it:stub}{cmd:_eM}, {it:stub}{cmd:_eY}{p_end}
{synoptline}
{p 4 6 2}
* {cmd:mediator()}, {cmd:treatment()}, and {cmd:instrument()} are required.{p_end}
{p 4 6 2}
{it:indepvars} may contain factor variables ({help fvvarlist}) and time-series
operators ({help tsvarlist}). Factor-variable ({cmd:i.}) terms enter as
drop-base dummies. Exactly collinear controls are dropped (in order, with a
note), and the controls are partialed out by the Frisch-Waugh-Lovell theorem and
not reported. Wildcards and ranges in {it:indepvars} (for example {cmd:x*} or
{cmd:x1-x8}) never pick up the outcome, treatment, mediator, instruments, or the
{cmd:absorb()}, cluster, weight, and coordinate variables; those matches are
excluded with a note. Naming the outcome, treatment, mediator, or an instrument
among the controls is an error.{p_end}
{p 4 6 2}
{cmd:aweight}s and {cmd:pweight}s are allowed; see {help weight} and the
{help ivmediate##remarks:Remarks} below.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:ivmediate} decomposes the effect of a single treatment {it:T} on an
outcome {it:Y} into a {bf:direct} effect and an {bf:indirect} (mediation)
effect operating through a mediator {it:M}, using one or more excluded
instruments {it:Z} for {it:T} -- {it:without} requiring a separate instrument
for the mediator. Identification rests on the
{bf:Mediated Confounding Assumption} (MCA): the treatment- and
outcome-equation structural errors are uncorrelated. The framework and
identifying assumption are developed in section 2 (the IV mediation framework)
and section 3 (the Mediated Confounding Assumption); the linear identification
and 2SLS estimator are in section 4.

{pstd}
The command reports the standard IV total effect {it:tau_total} (the coefficient
on {it:T} in the 2SLS regression of {it:Y} on {it:T}), the direct effect
{it:tau}, the mediator response {it:gamma}, the outcome-on-mediator coefficient
{it:theta}, and the indirect effect {it:theta*gamma}, each with a Wald interval.
With one instrument the identity {it:tau} + {it:theta*gamma} = {it:tau_total}
holds algebraically; the command checks it numerically and warns if the relative
error exceeds 1e-8. With more than one instrument the equations are estimated one
by one by overidentified 2SLS, and the identity holds only asymptotically (a note
is printed when the relative gap exceeds 1e-2).

{pstd}
Alongside the point estimates {cmd:ivmediate} reports the first stage of {it:T}
on {it:Z} (pi1), first-stage diagnostics (the Olea-Pflueger [2013] effective F
for the treatment and the conditional-mediator first stages, with the
simplified-test critical value at the estimated effective degrees of freedom
{it:K_eff}), Anderson-Rubin weak-instrument-robust confidence sets (by default),
a {bf:kappa}-based sensitivity analysis (breakeven confounding, an inference
breakeven, and, with {opt kappatable}, a table of kappa-adjusted estimates), and,
with {opt r2sensitivity}, an {bf:R-squared}-based sensitivity analysis
(robustness values RV and XRV in the sense of Cinelli and Hazlett [2020]); both
are always computed and stored in {cmd:e()}. Inference uses heteroskedasticity-robust (the default),
classical, cluster-robust, or Conley spatial standard errors, optionally
complemented by a bootstrap; a post-double-selection LASSO path selects controls
from a basis expansion while keeping the diagnostics.


{marker papermap}{...}
{title:Where to find the theory}

{pstd}
Each of the main features is developed in the companion paper:

{synoptset 37}{...}
{synopt:{it:Feature}}{it:Paper location}{p_end}
{synoptline}
{synopt:Decomposition and identification}section 4 (Identification and Estimation under Linearity), esp. 4.1-4.2{p_end}
{synopt:Mediated Confounding Assumption}section 3 (The Mediated Confounding Assumption){p_end}
{synopt:Olea-Pflueger effective F}section 4.2, Remark (Conditional First-Stage Strength); section 7.3, Experiment 2{p_end}
{synopt:Anderson-Rubin confidence sets}section 4.2 (weak-instrument-robust inference){p_end}
{synopt:kappa-based sensitivity}section 5.1 (Sensitivity for the Mediator Effect); Appendix F.1, F.4{p_end}
{synopt:R-squared sensitivity (RV, XRV)}section 5.1; Appendix F.2-F.3 (R-squared-Based Sensitivity){p_end}
{synopt:Post-double-selection LASSO}section 4.2, Remark (High-dimensional controls); Appendix H{p_end}
{synoptline}
{p 4 6 2}
Section and appendix numbers refer to the working-paper draft distributed with
the replication materials; titles are given so the cross-references survive
minor renumbering.{p_end}


{marker whatsnew}{...}
{title:What's new in 2.0}

{pstd}
Version 2.0.0 replaces {cmd:ivmediate} 1.0.4 (SSC, 2020; Dippel, Ferrara, and
Heblich 2020) under the same name. The syntax {it:depvar} [{it:indepvars}]{cmd:,}
{opt mediator()} {opt treatment()} {opt instrument()} is kept, but the estimator,
the defaults, and the output change; results obtained with version 1 should be
re-run. Relative to 1.0.4:

{p 8 10 2}o  {bf:Defaults:} heteroskedasticity-robust (HC1) standard errors and
Anderson-Rubin confidence sets (version 1 reported unadjusted standard errors);
{cmd:vce(unadjusted)} gives classical standard errors and {cmd:noar} skips the
sets.{p_end}
{p 8 10 2}o  {bf:{cmd:full} removed:} the option is no longer accepted.{p_end}
{p 8 10 2}o  {bf:Instruments:} {opt instrument()} takes {it:K} >= 1 instruments.
Instrument strength is summarized by the Olea-Pflueger effective F, with critical
values at the estimated {it:K_eff}, instead of the Kleibergen-Paap F.{p_end}
{p 8 10 2}o  {bf:Fixed effects:} {opt absorb()} takes a varlist of one or two
categorical variables (one- or two-way fixed effects).{p_end}
{p 8 10 2}o  {bf:Weights and standard errors:} {cmd:aweight}s and {cmd:pweight}s;
cluster-robust (Liang-Zeger CR1) and Conley spatial standard errors.{p_end}
{p 8 10 2}o  {bf:LASSO and bootstrap:} post-double-selection LASSO over a basis
expansion of the controls ({opt lasso}), and bootstrap standard errors and
confidence intervals ({opt bootstrap}).{p_end}
{p 8 10 2}o  {bf:Sensitivity analysis:} kappa-based (the breakeven kappa* and
kappa*_CI, always reported; the kappa-adjusted table with {opt kappatable}) and
R-squared-based (the robustness values RV and XRV, printed with
{opt r2sensitivity}) sensitivity to violations of the Mediated Confounding
Assumption.{p_end}
{p 8 10 2}o  {bf:No dependencies:} {cmd:ranktest} is no longer needed; the
computations run in Mata.{p_end}
{p 8 10 2}o  {bf:Estimation:} equation-by-equation 2SLS rather than a joint
one-step GMM system. With one instrument and {cmd:vce(unadjusted)} the total,
direct, and indirect point estimates coincide with version 1 to machine
precision, but the standard errors differ because version 1 reported GMM system
standard errors.{p_end}


{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{opth mediator(varname)}, {opth treatment(varname)} specify the single mediator
and single treatment. Both are required. The outcome, treatment, and mediator
must be three different variables.

{phang}
{opth instrument(varlist)} specifies one or more excluded instruments for the
treatment. Required. Instruments must differ from the outcome, treatment, and
mediator, and may not be listed twice. With one instrument the reported
{it:F(T,Z|X)} and {it:F(M,Z|T,X)} are the effective F (= t-squared); with more
than one, the joint Wald F is stored in {cmd:e(F_TZ)} and {cmd:e(F_MZ_T)}, the
Olea-Pflueger effective F in {cmd:e(F_TZ_op)} and {cmd:e(F_MZ_T_op)}, and the
per-instrument t-squared values of the conditional first stage in
{cmd:e(F_MZ_T_individual)}. See section 4 for identification and estimation under
overidentification.

{phang}
{opth absorb(varlist)} absorbs one or two categorical variables (numeric or
string) as fixed effects by within-group demeaning (Frisch-Waugh-Lovell). String
values that are empty or contain only blanks count as missing. Singleton groups
are dropped iteratively across dimensions before estimation; a sample in which
every group is a singleton is an error. More than two dimensions is an error.

{dlgtab:Standard errors}

{phang}
{cmd:vce(}{it:vcetype}{cmd:)} selects the variance estimator. {it:vcetype} is
{cmd:robust} (HC1; the default; also {cmd:r}, {cmd:rob}),
{cmd:unadjusted} (homoskedastic; also {cmd:un}, {cmd:ols}, {cmd:classical}),
{cmd:cluster} {it:clustvar} (Liang-Zeger CR1, one-way; also {cmd:cl}, {cmd:clu},
{cmd:clus}; the cluster variable may be numeric or string), or {cmd:conley}
(spatial HAC; also {cmd:spatial}). Other words, or further tokens inside
{cmd:vce()}, are an error. Under {cmd:pweight}s an {cmd:unadjusted} request is
upgraded to {cmd:robust}. Cluster SEs need at least 2 clusters.

{phang}
{opt latitude(varname)}, {opt longitude(varname)} give the coordinates used by
{cmd:vce(conley)}; both are required with it. Under {cmd:metric(haversine)} they
are latitude and longitude in decimal degrees (values outside [-90, 90] and
[-180, 180] are an error); under {cmd:metric(euclidean)} they are y- and
x-coordinates in the units of {opt cutoff()}. Specifying {opt latitude()},
{opt longitude()}, {opt cutoff()}, {opt kernel()}, or {opt metric()} without
{cmd:vce(conley)} is an error.

{phang}
{opt cutoff(#)} is the distance beyond which observations are treated as
spatially independent (kilometers under {cmd:haversine}, coordinate units under
{cmd:euclidean}). Required with {cmd:vce(conley)}. Conley SEs are sensitive to
this choice; report results at several cutoffs. With the uniform kernel a cutoff
that puts every pair of observations within range is an error (the Conley meat
is zero by construction), and a cutoff that covers more than half of all pairs
prints a warning. The pairs within the cutoff are stored sparsely; computing the
distances takes time of order {it:N}-squared, and more than 100 million pairs
within the cutoff is an error.

{phang}
{opt kernel(string)} sets the spatial kernel: {cmd:uniform} (equal weight
within {opt cutoff()}; the default) or {cmd:bartlett} (linear taper to zero at
{opt cutoff()}). Neither kernel matrix is guaranteed positive semidefinite in two
dimensions (the triangular radial function is positive definite only in one
dimension), so a Conley variance can be negative; see
{help ivmediate##remarks:Remarks}. {opt kernel()} and {opt metric()} are not
case sensitive.

{phang}
{opt metric(string)} sets the distance metric: {cmd:haversine} (great-circle
distance on a sphere of radius 6,371 km; the default) or {cmd:euclidean}
(planar).

{dlgtab:Weak-instrument diagnostics and robust sets}

{phang}
{opt optau(#)} and {opt opalpha(#)} set the bias fraction {it:tau} and size
{it:alpha} of the Olea-Pflueger weak-instrument test; both must lie strictly
between 0 and 1. The critical value is the simplified-test value at the estimated
effective degrees of freedom {it:K_eff} (Montiel Olea and Pflueger 2013;
{cmd:weakivtest}'s c_simp), computed for any {it:tau} and {it:alpha}; see
{help ivmediate##remarks:Remarks}. With one instrument {it:K_eff} = 1 and the
critical value at {it:tau} = 0.10, {it:alpha} = 0.05 is 23.109. The conditional
first stage and the effective-F diagnostic are motivated in section 4.2 (Remark,
Conditional First-Stage Strength) and calibrated in section 7.3 (Experiment 2).

{phang}
{opt noar} skips the Anderson-Rubin weak-instrument-robust confidence sets, which
are computed by default for the mediator effect {it:theta}, the
treatment-on-mediator effect {it:gamma}, and the total effect, at the confidence
{opt level()} and under the same SE regime as the Wald inference. {opt ar} is
accepted and has no effect. With one instrument each set is obtained by exact
Fieller-type test inversion; with several instruments the AR statistic is
computed in closed form, the accepted segments are located on a grid and their
bounds refined by Brent's method. A set may be a bounded interval, a pair of rays,
the whole line, a union of segments ({it:K}>1), empty, or unavailable (negative
Conley variance); the shape is stored in {cmd:e(ar_theta_shape)},
{cmd:e(ar_gamma_shape)}, and {cmd:e(ar_total_shape)}. An unbounded or
disconnected set is the expected outcome under weak identification -- the data
cannot bound the parameter at that confidence level; the output says so only when
a set is unbounded or disconnected. AR inference is recommended when the effective
F falls below its critical value (section 4.2); see Anderson and Rubin (1949).

{phang}
{opt arpoints(#)} (at least 3) and {opt arwidth(#)} (positive) control the AR
search grid when there is more than one instrument: the number of grid points and
the half-width of the grid in standard-error units around the point estimate. The
range {opt arwidth()} x SE is also the range on which a Conley variance is checked
for negativity, also with one instrument. With one instrument the set itself is
exact and does not depend on the grid.

{dlgtab:Bootstrap}

{phang}
{opt bootstrap} adds bootstrap standard errors and confidence intervals, printed
in a separate block after the analytic results, for the total, direct, mediator,
treatment-on-mediator, and indirect effects and for kappa* and rho*_TY (and,
with {opt bdiagnostics}, RV and XRV). The estimates table, {cmd:e(b)}, and
{cmd:e(V)} stay analytic. The bootstrap is not available with {cmd:lasso}
(bootstrapping the selection under-covers; use the analytic post-double-selection
inference).

{phang}
{opt breps(#)} sets the number of replications (at least 50; 999 recommended).
{opt bmethod(string)} chooses the resampling scheme: {cmd:block} (cluster
resampling; the default with {cmd:vce(cluster)} and only available then) or
{cmd:pairs} (resampling of observations; the default otherwise; valid for
independent observations only). {opt bcitype(string)} chooses the interval:
{cmd:basic} (the default), {cmd:percentile}, or {cmd:studentized}
(bootstrap-t with each replicate's analytic SE). kappa*, rho*_TY, RV, and XRV
have no analytic SE, so under {cmd:studentized} they get the basic interval.
{opt bmethod()} and {opt bcitype()} are not case sensitive.

{phang}
{opt bseed(#)} sets the random-number seed (a nonnegative integer). The command
saves the current random-number state ({cmd:c(rngstate)}), runs {cmd:set seed},
draws the replicates sequentially, and restores the saved state afterward, also
after an error: the user's random-number stream is unchanged, and two different
seeds give different replicates.

{phang}
{opt bootidx(matname)} supplies the resamples instead of drawing them, one row per
replicate; the number of rows sets the number of replications ({opt breps()} and
{opt bseed()} are then not used, and {cmd:e(boot_seed)} is missing). Entries are
1-based integers referring to the estimation sample, that is, the observations
with complete data and positive weight (before singleton groups are dropped), in
data order. Under {cmd:pairs} the matrix has {it:n} columns and each entry is an
observation of that sample. Under {cmd:block} the matrix has {it:G} columns and
each entry is a cluster, numbered in order of first appearance in that sample. A
matrix that does not exist, has the wrong number of columns, or contains entries
outside 1 to {it:n} (or {it:G}) is an error.

{phang}
{opt bdiagnostics} also bootstraps RV(theta) and XRV(theta).

{dlgtab:Post-double-selection LASSO}

{phang}
{opt lasso} selects controls for {it:Y}, {it:T}, {it:M}, and each instrument from
a basis expansion of the supplied controls by the quadruple-LASSO variant of the
Belloni, Chernozhukov, and Hansen (2014) post-double-selection procedure, then
estimates the structural equations on the union of the selected terms (and any
{opt lkeep()} terms). The mediator is selected on the expansion alone, not on the
expansion and {it:T} (Appendix H, Algorithm 1). The method is stated in section
4.2 (Remark, High-dimensional controls via post-double selection) and developed in
Appendix H. With {cmd:lasso}, {it:indepvars} may contain continuous variables and
{cmd:i.} factor terms, but no interactions (the expansion builds them).

{phang}
{opt lpenalty(string)} selects the penalty rule. Only {cmd:plugin} (the BCH
plug-in penalty; the default) is supported. A cross-validated penalty is not
offered: it depends on the random-number generator and cannot be reproduced
across platforms.

{phang}
{opt expandx(string)} lists the operators that build the expansion of the
continuous controls: any of {cmd:square}, {cmd:cube}, and {cmd:interact_with_X}
(pairwise products of the continuous controls), or {cmd:none}; default
{cmd:square cube}. Powers apply only to variables with more than two distinct
values. {opt expanddum(string)} lists the operators involving the factor dummies:
any of {cmd:interact_with_X} (dummy x continuous control) and
{cmd:interact_with_dummies} (dummy x dummy of different factors), or {cmd:none};
default {cmd:interact_with_X}. Operators are not case sensitive and are validated
whenever given. The expansion contains, in this order: the continuous controls,
the dummies, the squares and cubes, the products of continuous controls, the
dummy x control products, and the dummy x dummy products. Expansion columns that
are (nearly) constant on the sample are dropped before selection; collinear
columns are kept in the expansion, and exactly collinear {it:selected} controls
are dropped (with a note) before the post-selection 2SLS.

{phang}
{opt lkeep(patterns)} names expansion terms that are always kept: they are
partialed out of every LASSO and added to the selected set. Each pattern is
matched against the expansion term names with {help strmatch()}: {cmd:*} matches
any string and {cmd:?} one character. Names are as displayed: {cmd:turn},
{cmd:turn^2}, {cmd:turn#trunk}, {cmd:2.rep78}, {cmd:2.rep78#turn}; a {cmd:c.}
prefix is ignored, so {cmd:lkeep(c.turn^2)} is the same as {cmd:lkeep(turn^2)}.
Because the match is on the name, {cmd:c1*} keeps {cmd:c1}, {cmd:c1^2},
{cmd:c1#c2}, but also {cmd:c10} and its terms, and it does not keep {cmd:c0#c1} or
{cmd:2.g#c1}; list such terms separately. A pattern that matches no term is an
error. {opt lkeep()} is ignored without {opt lasso}.

{phang}
{opt showsel} adds the per-equation selections to the output: the terms each of
the LASSOs (for {it:Y}, {it:T}, {it:M}, and each instrument) selected, with its
penalty level and whether its loading iterations converged ("converged in {it:k}
iterations" or "did not converge after {it:k} iterations").

{phang}
{opt lcslack(#)} (the slack {it:c} > 0, default 1.1), {opt lalphap(#)} (the tail
probability in (0, 1), default 0.1/ln({it:n})), {opt lniter(#)} (default 15), and
{opt ltol(#)} (default 1e-6) are advanced controls of the plug-in penalty
lambda = 2 {it:c} sqrt({it:n}) invnormal(1 - {it:alpha_p}/(2{it:p})); see
{help ivmediate##remarks:Remarks}. {opt llambda(#)} fixes lambda at {it:#}
instead (for replicating a published penalty); the loadings are then computed once
from the centered dependent variable.

{dlgtab:Reporting}

{phang}
{opt level(#)} sets the confidence level; default {cmd:level(95)} or as set by
{helpb set level}. It governs the Wald, bootstrap, and AR intervals and the RV
threshold. Labels show the level as typed, e.g., 97.5%.

{phang}
{opt nostars} suppresses the significance stars; {opt starcutoffs(# # #)} sets
the three two-sided p-value thresholds (in any order; the strictest is assigned
three stars). The legend reads "*** p<0.01, ** p<0.05, * p<0.1" by default.

{phang}
{opt kappatable} prints the table of kappa-adjusted estimates below the kappa
summary lines: theta adjusted for kappa0 = 0, 10%, 25%, 50%, 75%, 90%, and 100% of
kappa*, with Wald intervals. When the effective F of the conditional first stage is
below its OP critical value ({cmd:e(F_MZ_T_above_crit)} = 0), a note below the
table says that c_TM is estimated imprecisely and, when the Anderson-Rubin sets are
shown, that the AR set for theta is robust to weak instruments. The table is always
stored in {cmd:e(kappa_adj)}.

{phang}
{opt r2sensitivity} prints the R-squared-based sensitivity section: the reduced-form
F (= t-squared), RV(theta), XRV(theta), and their reduced-form and first-stage
components. The values are always stored in {cmd:e()}.

{phang}
{opt residuals(stub)} saves, on the estimation sample, the residuals that feed the
kappa block: {it:stub}{cmd:_eT} (first stage of {it:T} on {it:Z}, after
partialing out the controls), {it:stub}{cmd:_eM} ({it:M} - gamma {it:T}, same
scale), and {it:stub}{cmd:_eY} (the conditional 2SLS residual, unweighted scale).


{marker remarks}{...}
{title:Remarks}

{pstd}
{bf:Estimation sample.} The sample is the set of observations with nonmissing
values on all model variables: the outcome, controls, mediator, treatment,
instruments, and -- where applicable -- the cluster variable, the {opt absorb()}
variables, the coordinates, and the weight. Empty and blank-only string values of
the cluster and {opt absorb()} variables count as missing. Observations with zero
weight are dropped next (with a note), then singleton FE groups are dropped
iteratively. {cmd:e(sample)} marks the sample used. Negative weights are an
error, and so are all weights zero.

{pstd}
{bf:Weights.} Weights enter by sqrt(w)-premultiplication; the results do not
depend on the scale of the weights. With {opt absorb()} the fixed effects are
removed with weighted group means, so absorbing a factor gives the same results as
including its dummies. The kappa block uses Hajek-form weighted moments
(Appendix F.4).

{pstd}
{bf:Fixed effects.} One-way FE are removed exactly. With two-way FE the command
computes the exact projection when the dimension with fewer levels has at most
3,000 levels (it eliminates the larger dimension and solves the remaining
Laplacian system, one level fixed per connected component); otherwise it uses
alternating projections with Irons-Tuck acceleration, which stop when the largest
change in one sweep is below 1e-12 times the standard deviation of the column
(at most 10,000 iterations; a warning is printed only if they do not converge).
For the degrees of freedom see {help ivmediate##conventions:Conventions}.

{pstd}
{bf:Controls, identification, and numerics.} Every regression is solved by a QR
decomposition of the design with columns scaled to unit length; for 2SLS the
fitted values come from a QR of the instrument matrix. A control is dropped as
collinear when its residual on the intercept and the earlier controls (after the
fixed effects) is below 1e-9 of its variation; the earlier member of a collinear
set is kept, as with {cmd:_rmcoll}. The command stops with an error when, after
the controls and fixed effects, an instrument has no variation, the instruments
are collinear, the treatment or the mediator has no variation, or the mediator is
collinear with the treatment, and when the number of observations does not exceed
the number of coefficients of the largest regression plus the absorbed degrees of
freedom.

{pstd}
{bf:Undefined variances.} The Conley sandwich is not guaranteed positive
semidefinite, and a variance can come out negative. Such a standard error is
reported as missing ("SE unavailable" in the table), the F statistics and
{it:K_eff} that use a variance block that is not positive definite are missing, an
Anderson-Rubin set whose variance turns negative on its search range is reported
as unavailable, and one warning is printed. If the variance of the total, direct,
or indirect effect is unavailable, {cmd:e(b)} is posted without {cmd:e(V)} and
{cmd:e(V_unavailable)} is 1; all other results are kept.

{pstd}
{bf:Olea-Pflueger critical values.} Let W2 = Q^(1/2) V Q^(1/2), with V the
variance of the instrument coefficients under the chosen SE regime and Q =
Z'Z/{it:n} for the instruments residualized on the controls. Then {it:K_eff} =
[tr W2]^2 (1 + 2x) / (tr(W2'W2) + 2x tr(W2) lambda_max(W2)), x = 1/{it:tau}, and
the critical value is {cmd:invnchi2(}{it:K_eff}{cmd:, }x{it:K_eff}{cmd:, }1 -
{it:alpha}{cmd:)}/{it:K_eff}. For {it:K_eff} = {it:K} this reproduces Table 1 of
Montiel Olea and Pflueger (2013); {cmd:F_TZ} and {cmd:F_MZ_T} each get their own
{it:K_eff}.

{pstd}
{bf:kappa-based sensitivity.} kappa* = theta c_TM is the breakeven value of
Cov(eps_T, eps_Y); kappa*_CI = (theta - sign(theta) z SE(theta)) c_TM is the
value at which the kappa-adjusted CI first includes 0. When the CI already
contains 0, kappa*_CI is flagged as uninformative
({cmd:e(kappa_star_CI_informative)} = 0). {cmd:e(F_MZ_T_above_crit)} is 1 when the
effective F of the conditional first stage is at least its OP critical value at
{opt optau()} and 0 otherwise; in the latter case c_TM is estimated imprecisely,
and with {opt kappatable} a note below the table says so (paper App. F.4, step 3).
When |rho*_TY| > 1, no feasible error correlation overturns theta.

{pstd}
{bf:R-squared-based sensitivity.} F (reduced form) and F (conditional first
stage) are the squared t statistics of the instrument in the regressions of {it:Y}
and {it:M} on [1, Z] after partialing out [1, X, T] and the fixed effects, under
the chosen SE regime ({cmd:e(F_reduced_form)}, {cmd:e(F_first_stage_cond)}). The
section printed by {opt r2sensitivity} omits the conditional first-stage F, which
repeats F(M, Z|T,X) of the first-stage diagnostics. RV and XRV use the residual
degrees of freedom of the regression of {it:Y} on [1, X, T, Z] (as
{cmd:sensemakr}). With more than one instrument they use the first instrument; the
output says which.

{pstd}
{bf:LASSO.} Each LASSO centers the dependent variable and the candidate columns at
their weighted means (an unpenalized intercept) after partialing out the
{opt lkeep()} terms, and uses the weights normalized to mean 1. The penalty is
lambda = 2 {it:c} sqrt({it:n}) invnormal(1 - {it:alpha_p}/(2{it:p})), with {it:n}
the number of observations (the Kish size (sum w)^2/sum w^2 with unequal weights)
also under clustering or Conley SEs, where the dependence enters through the
loadings (Belloni, Chernozhukov, Hansen, and Kozbur 2016; lassopack). The
default {it:alpha_p} = 0.1/ln({it:n}) uses the number of observations also under
clustering, where lassopack uses 0.1/ln({it:G}). The loadings follow the SE
regime (heteroskedastic, cluster, or Conley scores w x W x u-hat) and are updated
from the residuals of a post-LASSO regression on the selected terms, the
{opt lkeep()} terms, and an intercept, until the coefficients change by less
than {opt ltol()} or {opt lniter()} iterations have run. Each penalized fit is a
coordinate-descent port of scikit-learn's Lasso (tolerance 1e-6, at most 20,000
sweeps).

{pstd}
{bf:Bootstrap.} The bootstrap resamples the estimation sample (complete cases with
positive weight); block resampling draws its clusters, each drawn copy being a
separate cluster, and FE nested in the clusters get new levels for every copy (as
{cmd:bootstrap, cluster() idcluster()}), so the replicates keep the nesting.
Each replicate repeats the singleton pass, demeaning, and collinearity checks; a
replicate that fails is counted ({cmd:e(boot_n_failed)}; a warning is printed
above 10%). For {cmd:basic} and {cmd:percentile} intervals without
{opt bdiagnostics} the replicates compute point estimates only (fast, also with
Conley SEs); {cmd:studentized} intervals and {opt bdiagnostics} refit fully. Under
{cmd:vce(conley)} the resampling does not reproduce continuous spatial
dependence; consider clustering on spatial cells instead. The bootstrap block also
reports kappa*_CI computed with the bootstrap SE of theta.

{pstd}
{bf:Cross-language validation.} The Mata engine ports the Python reference
implementation ({cmd:ivmediate} for Python, version 2.0.0). The harness in
{cmd:validation/} of the
{browse "https://github.com/aferrara87/ivmediate":GitHub repository} (run by
{cmd:validation/stata/run_stata.do}) compares the two on
the configurations of {cmd:validation/configs.json}: classical, robust, cluster,
and Conley SEs; analytic and probability weights; one- and two-way FE (nested,
disconnected, weighted); dummies; K = 1 and 2; the LASSO variants; the bootstrap
with shared resamples; and inputs that must stop with an error. Stata agrees with
Python on every configuration it can express (relative tolerance 1e-8); two are
not expressible in Stata (an infinite value, and per-variable expansion rules).


{marker conventions}{...}
{title:Conventions}

{p 4 6 2}
o  Wald confidence intervals use normal critical values z(1 - alpha/2);
Anderson-Rubin sets use chi2({it:K}) critical values.{p_end}
{p 4 6 2}
o  Absorbed fixed effects count as the rank of the FE dummies beyond the
intercept: {it:G} - 1 one-way; {it:G1} + {it:G2} - {it:c} - 1 two-way, with {it:c}
the number of connected components of the bipartite FE graph ({cmd:e(df_fe)}).
Classical and HC1 SEs count all of them. With {cmd:vce(cluster)}, FE dimensions
nested in the clusters (every FE group lies within one cluster) are not counted in
the CR1 factor (reghdfe/fixest convention); when one of two dimensions is nested,
the other counts {it:G} - 1 ({cmd:e(df_fe_vce)}, {cmd:e(fe_nested)}).{p_end}
{p 4 6 2}
o  CR1: {it:G}/({it:G}-1) x ({it:n}-1)/({it:n}-{it:k}-df); HC1:
{it:n}/({it:n}-{it:k}-df); Conley: ({it:n}-1)/({it:n}-{it:k}-df); classical:
sigma-squared with {it:n}-{it:k}-df.{p_end}
{p 4 6 2}
o  Olea-Pflueger critical values: the simplified test at the estimated {it:K_eff}
({it:K_eff} = 1 for one instrument).{p_end}
{p 4 6 2}
o  RV and XRV with {it:K} >= 2 instruments use the first instrument.{p_end}


{marker examples}{...}
{title:Examples}

{pstd}Setup (illustrative roles on the auto data; not a causal claim){p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}

{pstd}Single-instrument mediation analysis (robust SEs and Anderson-Rubin sets by
default){p_end}
{phang2}{cmd:. ivmediate price foreign turn, treatment(weight) mediator(length) instrument(displacement)}{p_end}

{pstd}Cluster-robust standard errors{p_end}
{phang2}{cmd:. ivmediate price foreign turn, treatment(weight) mediator(length) instrument(displacement) vce(cluster rep78)}{p_end}

{pstd}Absorbed fixed effects and analytic weights{p_end}
{phang2}{cmd:. ivmediate price turn [aw=turn], treatment(weight) mediator(length) instrument(displacement) absorb(rep78)}{p_end}

{pstd}Two instruments: critical values at the estimated K_eff and AR sets from the
exact statistic (on these data the conditional first stage is weak, so the theta
set is unbounded -- the diagnostic working as intended){p_end}
{phang2}{cmd:. ivmediate price foreign turn, treatment(weight) mediator(length) instrument(displacement gear_ratio)}{p_end}

{pstd}Conley spatial standard errors (illustrative coordinates){p_end}
{phang2}{cmd:. ivmediate price foreign, treatment(weight) mediator(length) instrument(displacement) vce(conley) latitude(turn) longitude(trunk) cutoff(20) kernel(bartlett) metric(euclidean)}{p_end}

{pstd}Bootstrap confidence intervals in addition to the analytic ones{p_end}
{phang2}{cmd:. ivmediate price foreign turn, treatment(weight) mediator(length) instrument(displacement) bootstrap breps(999) bcitype(percentile)}{p_end}

{pstd}Post-double-selection LASSO over a basis expansion, forcing one term and
showing the per-equation selections{p_end}
{phang2}{cmd:. ivmediate price c.turn c.trunk i.rep78, treatment(weight) mediator(length) instrument(displacement) lasso lkeep(turn^2) showsel}{p_end}

{pstd}The kappa-adjusted table and the R-squared-based sensitivity section{p_end}
{phang2}{cmd:. ivmediate price foreign turn, treatment(weight) mediator(length) instrument(displacement) kappatable r2sensitivity}{p_end}

{pstd}Classical SEs and 90% intervals without the Anderson-Rubin sets{p_end}
{phang2}{cmd:. ivmediate price foreign turn, treatment(weight) mediator(length) instrument(displacement) vce(unadjusted) level(90) noar}{p_end}

{pstd}Results{p_end}
{phang2}{cmd:. ereturn list}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:ivmediate} stores the following in {cmd:e()}:

{synoptset 30 tabbed}{...}
{syntab:Scalars}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(V_unavailable)}}1 if a variance in {cmd:e(V)} is unavailable (then {cmd:e(V)} is not posted){p_end}
{synopt:{cmd:e(k_inst)}}number of instruments {it:K}{p_end}
{synopt:{cmd:e(k_cov)}}number of controls used (after collinear ones are dropped; after selection with {cmd:lasso}){p_end}
{synopt:{cmd:e(k_dummies)}}number of factor dummies among the controls (before selection){p_end}
{synopt:{cmd:e(df_fe)}}absorbed FE degrees of freedom (rank of the FE dummies){p_end}
{synopt:{cmd:e(df_fe_vce)}}absorbed FE degrees of freedom counted in the variance{p_end}
{synopt:{cmd:e(N_fegroups)}}total number of FE groups (all dimensions){p_end}
{synopt:{cmd:e(N_singleton)}}singleton observations dropped{p_end}
{synopt:{cmd:e(N_clust)}}number of clusters {it:G} (with {cmd:vce(cluster)}){p_end}
{synopt:{cmd:e(spatial_cutoff)}}spatial cutoff (with {cmd:vce(conley)}){p_end}
{synopt:{cmd:e(spatial_share)}}share of pairs within the cutoff (with {cmd:vce(conley)}){p_end}
{synopt:{cmd:e(spatial_negvar)}}1 if a reported Conley variance was negative{p_end}
{synopt:{cmd:e(total)}, {cmd:e(se_total)}}total effect and its SE{p_end}
{synopt:{cmd:e(tau)}, {cmd:e(se_tau)}}direct effect and its SE{p_end}
{synopt:{cmd:e(theta)}, {cmd:e(se_theta)}}outcome-on-mediator coefficient and SE{p_end}
{synopt:{cmd:e(gamma)}, {cmd:e(se_gamma)}}treatment-on-mediator coefficient and SE{p_end}
{synopt:{cmd:e(indirect)}, {cmd:e(se_indirect)}}indirect effect and its delta-method SE{p_end}
{synopt:{cmd:e(cov_gamma_theta)}}cross-equation covariance of gamma and theta{p_end}
{synopt:{cmd:e(z_crit)}}normal critical value used for the Wald CIs{p_end}
{synopt:{cmd:e(decomp_check)}}tau + theta*gamma{p_end}
{synopt:{cmd:e(decomp_error)}}|tau + theta*gamma - total|{p_end}
{synopt:{cmd:e(decomp_error_rel)}}decomp_error / max(|total|, |tau| + |theta*gamma|){p_end}
{synopt:{cmd:e(mepct)}}indirect as a percentage of the total effect (if the total is nonzero){p_end}
{synopt:{cmd:e(F_TZ)}, {cmd:e(F_TZ_op)}}treatment first-stage F (joint Wald; OP effective){p_end}
{synopt:{cmd:e(K_eff_TZ)}, {cmd:e(F_TZ_crit)}}its effective degrees of freedom and OP critical value{p_end}
{synopt:{cmd:e(F_MZ_T)}, {cmd:e(F_MZ_T_op)}}conditional first-stage F (joint Wald; OP effective){p_end}
{synopt:{cmd:e(K_eff_MZ_T)}, {cmd:e(F_MZ_T_crit)}}its effective degrees of freedom and OP critical value{p_end}
{synopt:{cmd:e(op_tau)}, {cmd:e(op_alpha)}}OP settings{p_end}
{synopt:{cmd:e(c_TM)}}Cov(eps_T, eps_M){p_end}
{synopt:{cmd:e(sigma_eT)}, {cmd:e(sigma_eY)}}residual SDs feeding rho*_TY{p_end}
{synopt:{cmd:e(kappa_star)}}breakeven confounding kappa*{p_end}
{synopt:{cmd:e(kappa_star_CI)}}inference breakeven kappa*_CI{p_end}
{synopt:{cmd:e(kappa_star_CI_informative)}}1 if the CI of theta excludes 0, else 0 (missing if SE unavailable){p_end}
{synopt:{cmd:e(rho_TY_star)}}breakeven expressed as a correlation{p_end}
{synopt:{cmd:e(F_MZ_T_above_crit)}}1 if the effective F of the conditional first stage is at least its OP critical value, else 0 (missing if either is unavailable){p_end}
{synopt:{cmd:e(F_reduced_form)}, {cmd:e(F_first_stage_cond)}}F (= t-squared) of the reduced form and the conditional first stage{p_end}
{synopt:{cmd:e(f_reduced_form)}, {cmd:e(f_first_stage_cond)}}the same (deprecated names){p_end}
{synopt:{cmd:e(RV_theta)}, {cmd:e(XRV_theta)}}robustness values for theta{p_end}
{synopt:{cmd:e(rv_rf)}, {cmd:e(rv_fs)}, {cmd:e(xrv_rf)}, {cmd:e(xrv_fs)}}their reduced-form and first-stage components{p_end}
{synopt:{cmd:e(level)}}confidence level{p_end}
{synopt:{cmd:e(ar)}}1 if the AR sets were computed{p_end}
{synopt:{cmd:e(lasso)}}1 if the LASSO path was used{p_end}
{synopt:{cmd:e(bootstrap)}}1 if bootstrapped{p_end}

{pstd}With the Anderson-Rubin sets ({it:x} is {cmd:theta}, {cmd:gamma}, or {cmd:total}):{p_end}
{synopt:{cmd:e(ar_level)}}confidence level of the sets{p_end}
{synopt:{cmd:e(ar_df)}, {cmd:e(ar_cv)}}degrees of freedom and chi2 critical value{p_end}
{synopt:{cmd:e(ar_grid_points)}, {cmd:e(ar_grid_halfwidth)}}search-grid settings{p_end}
{synopt:{cmd:e(ar_}{it:x}{cmd:_stat0)}, {cmd:e(ar_}{it:x}{cmd:_p)}}AR statistic and p-value at the null of zero (missing if the variance at 0 is not positive definite){p_end}
{synopt:{cmd:e(ar_}{it:x}{cmd:_contains0)}}1 if the set contains 0 (missing if unavailable){p_end}
{synopt:{cmd:e(ar_}{it:x}{cmd:_lower)}, {cmd:e(ar_}{it:x}{cmd:_upper)}}interval endpoints (shape {cmd:interval}){p_end}
{synopt:{cmd:e(ar_}{it:x}{cmd:_left_upper)}, {cmd:e(ar_}{it:x}{cmd:_right_lower)}}inner edges of the rays{p_end}

{pstd}With {cmd:bootstrap}:{p_end}
{synopt:{cmd:e(boot_B)}}replications{p_end}
{synopt:{cmd:e(boot_n_failed)}}replications that failed{p_end}
{synopt:{cmd:e(boot_n_est)}}size of the resampled estimation sample{p_end}
{synopt:{cmd:e(boot_seed)}}random-number seed (missing with {opt bootidx()}){p_end}
{synopt:{cmd:e(boot_alpha)}}1 - level/100{p_end}
{synopt:{cmd:e(boot_se_total)} ... {cmd:e(boot_se_indirect)}}bootstrap SEs of the five effects{p_end}

{pstd}With {cmd:lasso}:{p_end}
{synopt:{cmd:e(lasso_n_total)}}expansion terms offered to selection{p_end}
{synopt:{cmd:e(lasso_n_selected)}}terms in the union of the selections (incl. {opt lkeep()} terms){p_end}
{synopt:{cmd:e(lasso_n_keep)}}terms forced by {opt lkeep()}{p_end}

{synoptset 30 tabbed}{...}
{syntab:Macros}
{synopt:{cmd:e(cmd)}}{cmd:ivmediate}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(version)}}{cmd:2.0.0}{p_end}
{synopt:{cmd:e(title)}}{cmd:Linear IV Mediation Analysis}{p_end}
{synopt:{cmd:e(depvar)}}name of the outcome{p_end}
{synopt:{cmd:e(treat)}}name of the treatment{p_end}
{synopt:{cmd:e(med)}}name of the mediator{p_end}
{synopt:{cmd:e(inst)}}name(s) of the instrument(s){p_end}
{synopt:{cmd:e(controls)}}controls used (after collinear ones are dropped){p_end}
{synopt:{cmd:e(dropped_collinear)}}controls dropped as collinear (if any){p_end}
{synopt:{cmd:e(vce)}}{cmd:unadjusted}, {cmd:robust}, {cmd:cluster}, or {cmd:conley}{p_end}
{synopt:{cmd:e(vcetype)}}empty, {cmd:Robust}, or {cmd:Conley} (title used by {cmd:ereturn display}){p_end}
{synopt:{cmd:e(se_type)}}{cmd:classical}, {cmd:robust}, {cmd:cluster}, or {cmd:spatial} (the names used by the other implementations){p_end}
{synopt:{cmd:e(clustvar)}}cluster variable (with {cmd:vce(cluster)}){p_end}
{synopt:{cmd:e(wtype)}, {cmd:e(wexp)}}weight type and expression (if weighted){p_end}
{synopt:{cmd:e(absorb)}}absorbed FE variable(s) (if any){p_end}
{synopt:{cmd:e(latitude)}, {cmd:e(longitude)}}coordinate variables (with {cmd:vce(conley)}){p_end}
{synopt:{cmd:e(spatial_kernel)}, {cmd:e(spatial_metric)}}spatial kernel and metric (with {cmd:vce(conley)}){p_end}
{synopt:{cmd:e(ar_}{it:x}{cmd:_shape)}}{cmd:interval}, {cmd:rays}, {cmd:rays_left_only}, {cmd:rays_right_only}, {cmd:line}, {cmd:grid} ({it:K}>1), {cmd:empty}, or {cmd:unavailable}{p_end}
{synopt:{cmd:e(ar_}{it:x}{cmd:_method)}}{cmd:fieller} ({it:K}=1) or {cmd:grid}{p_end}
{synopt:{cmd:e(boot_method)}, {cmd:e(boot_citype)}}bootstrap method and CI type (if bootstrapped){p_end}
{synopt:{cmd:e(lasso_penalty)}}LASSO penalty rule (with {cmd:lasso}){p_end}
{synopt:{cmd:e(lasso_expandx)}, {cmd:e(lasso_expanddum)}}expansion operators used{p_end}
{synopt:{cmd:e(lasso_names)}}names of all expansion terms, in expansion order{p_end}
{synopt:{cmd:e(lasso_selected)}}names of the selected terms{p_end}
{synopt:{cmd:e(lasso_methods)}}per LASSO: {cmd:plugin}, {cmd:override} ({opt llambda()}), {cmd:keep-only}, or {cmd:trivial}{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}, or {cmd:b} if {cmd:e(V_unavailable)} = 1{p_end}

{synoptset 30 tabbed}{...}
{syntab:Matrices}
{synopt:{cmd:e(b)}}(total_effect, direct_effect, indirect_effect){p_end}
{synopt:{cmd:e(V)}}diagonal matrix of the three effect variances (not posted if one is unavailable){p_end}
{synopt:{cmd:e(pi1)}, {cmd:e(se_pi1)}}first stage of {it:T} on {it:Z} given the controls: coefficients and SEs (1 x {it:K}){p_end}
{synopt:{cmd:e(F_MZ_T_individual)}}per-instrument t-squared of the conditional first stage ({it:K} x 1){p_end}
{synopt:{cmd:e(kappa_adj)}}kappa-adjusted estimates (printed with {opt kappatable}): kappa0, fraction of kappa*, theta_adj, CI bounds{p_end}
{synopt:{cmd:e(fe_nested)}}1 per FE dimension nested in the clusters (with {cmd:vce(cluster)} and {opt absorb()}){p_end}
{synopt:{cmd:e(ar_}{it:x}{cmd:_seg)}}AR segments ({it:K}>1), one row per segment: {cmd:lower}, {cmd:upper}, {cmd:open_left}, {cmd:open_right} (1 if the segment reaches the edge of the search grid){p_end}
{synopt:{cmd:e(boot_point)}, {cmd:e(boot_se)}, {cmd:e(boot_ci)}}bootstrap point, SE, and CI of each bootstrapped quantity{p_end}
{synopt:{cmd:e(boot_analytic_se)}}analytic SEs of the same quantities (missing where none exists){p_end}
{synopt:{cmd:e(boot_dist)}}the replicates ({it:B} x quantities; missing rows are failed replicates){p_end}
{synopt:{cmd:e(lasso_sel_idx)}}positions of the selected terms in {cmd:e(lasso_names)}{p_end}
{synopt:{cmd:e(lasso_S)}, {cmd:e(lasso_B)}}per-LASSO selection indicators and coefficients ({it:p} x (3+{it:K}): Y, T, M, Z1, ...){p_end}
{synopt:{cmd:e(lasso_lambda)}, {cmd:e(lasso_n_iter)}, {cmd:e(lasso_converged)}}per-LASSO penalty, iterations, and convergence (1 x (3+{it:K})){p_end}

{synoptset 30 tabbed}{...}
{syntab:Functions}
{synopt:{cmd:e(sample)}}marks the estimation sample{p_end}

{pstd}
Note: {cmd:e(V)} is diagonal and holds the variances of the total, direct, and
indirect effects that correspond to the reported standard errors; the three
estimates come from different equations, so {cmd:e(V)} is not a joint covariance
matrix. The indirect-effect variance is the delta-method variance of
{it:theta*gamma}, using the cross-equation covariance. The bootstrapped quantities
are, in order, total, tau, theta, gamma, indirect, kappa*, rho*_TY, and with
{opt bdiagnostics} RV and XRV.


{marker citation}{...}
{title:Citation}

{pstd}
If you use {cmd:ivmediate}, please cite

{phang2}
Ferrara, Andreas, Robert Gold, Stephan Heblich, and Rodrigo Pinto (2026),
"Identifying Causal Mediation Effects with a Single Instrument," NBER Working
Paper No. [number], doi:[DOI].

{pstd}
and, for the original command, Dippel, Ferrara, and Heblich ({it:Stata Journal},
2020).


{marker references}{...}
{title:References}

{phang}
Anderson, T. W., and H. Rubin. 1949. Estimation of the parameters of a single
equation in a complete system of stochastic equations.
{it:Annals of Mathematical Statistics} 20: 46-63.

{phang}
Belloni, A., V. Chernozhukov, and C. Hansen. 2014. Inference on treatment
effects after selection among high-dimensional controls.
{it:Review of Economic Studies} 81: 608-650.

{phang}
Belloni, A., V. Chernozhukov, C. Hansen, and D. Kozbur. 2016. Inference in
high-dimensional panel models with an application to gun control.
{it:Journal of Business & Economic Statistics} 34: 590-605.

{phang}
Cinelli, C., and C. Hazlett. 2020. Making sense of sensitivity: Extending
omitted variable bias. {it:Journal of the Royal Statistical Society, Series B}
82: 39-67.

{phang}
Conley, T. G. 1999. GMM estimation with cross sectional dependence.
{it:Journal of Econometrics} 92: 1-45.

{phang}
Dippel, C., A. Ferrara, and S. Heblich. 2020. Causal mediation analysis in
instrumental-variables regressions. {it:Stata Journal} 20: 613-626. Software
update st0611_1, {it:Stata Journal} 21(1).

{phang}
Ferrara, A., R. Gold, S. Heblich, and R. Pinto. 2026. Identifying causal
mediation effects with a single instrument. NBER Working Paper No. [number],
doi:[DOI].

{phang}
Montiel Olea, J. L., and C. Pflueger. 2013. A robust test for weak instruments.
{it:Journal of Business & Economic Statistics} 31: 358-369.


{marker authors}{...}
{title:Authors}

{pstd}
Andreas Ferrara (maintainer), University of Pittsburgh and NBER,
{browse "mailto:a.ferrara@pitt.edu":a.ferrara@pitt.edu}{break}
Robert Gold, IfW Kiel and CESifo{break}
Stephan Heblich, University of Toronto and NBER{break}
Rodrigo Pinto, University of South Florida

{pstd}
Version 1 ({cmd:ivmediate} 1.0.4, 2020) by Christian Dippel, Andreas Ferrara, and
Stephan Heblich.

{pstd}
Questions and bug reports: the maintainer, or
{browse "https://github.com/aferrara87/ivmediate/issues":github.com/aferrara87/ivmediate/issues}.


{marker license}{...}
{title:License}

{pstd}
MIT license. Copyright (c) 2026 Andreas Ferrara, Robert Gold, Stephan Heblich,
Rodrigo Pinto. The license text is in the file LICENSE of the
{browse "https://github.com/aferrara87/ivmediate":GitHub repository}.


{marker also}{...}
{title:Also see}

{p 4 14 2}
Help:  {helpb ivregress}, {helpb areg}{p_end}
