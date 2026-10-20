function results = ivmediate(data, Y, T, M, Z, varargin)
%IVMEDIATE  IV mediation analysis with a single instrument (version 2.0.0).
%
%   results = ivmediate(data, Y, T, M, Z)
%   results = ivmediate(data, Y, T, M, Z, Name, Value, ...)
%
%   The IV mediation estimator of Ferrara, Gold, Heblich and Pinto,
%   "Identifying Causal Mediation Effects with a Single Instrument":
%   estimates, diagnostics, sensitivity analysis, weak-IV-robust
%   Anderson-Rubin sets, classical/robust/cluster/spatial standard errors,
%   weights, fixed effects, dummy expansion and post-double-selection LASSO;
%   bootstrap inference with ivmediate_bootstrap. Base MATLAB only (no
%   toolboxes); the Python package ivmediate is the reference
%   implementation.
%
%   Example:
%       data = ivmediate_read_csv('mydata.csv');      % or a table / struct
%       res = ivmediate(data, 'y', 'tvar', 'mvar', {'z1'}, ...
%                 'X', {'x1', 'x2'}, 'se', 'cluster', 'clusterid', 'clu');
%       ivmediate_print(res, 'title', 'My application');
%       bt = ivmediate_bootstrap(data, 'y', 'tvar', 'mvar', {'z1'}, ...
%                 'X', {'x1', 'x2'}, 'se', 'cluster', 'clusterid', 'clu');
%       ivmediate_print(res, 'boot', bt);  ivmediate_print_bootstrap(bt);
%   examples/becker_woessmann/bw.m in the repository reproduces the paper's
%   Table 6, column 1.
%
%   DATA
%   ----
%   A MATLAB table, a scalar struct of equal-length column vectors, or the
%   output of ivmediate_read_csv. Numeric and logical columns are numbers;
%   cellstr, string and categorical columns are text. Missing values: NaN in
%   numeric columns; '', whitespace-only, <missing> and <undefined> in text
%   columns. Only the columns an analysis uses are converted, so other
%   columns may have any type (datetime, duration, ...).
%
%   INPUTS AND OPTIONS
%   ------------------
%   Y, T, M: column names (text). Z: a name or a list of instrument names
%   (cellstr or string array; distinct, and different from Y, T, M).
%   Option names are case-insensitive. Defaults in brackets.
%     X                []   controls: cellstr / string array, or one
%                           whitespace-separated string; globs * ? [..]
%                           (case-sensitive) expand against the columns.
%                           Globs never pick up Y, T, M, Z or the fe /
%                           cluster / weight / coordinate columns; naming
%                           Y, T, M or Z explicitly is an error.
%     dummies          []   categorical variables expanded to indicators
%                           (all levels but the first sorted one); same
%                           input styles as X. Names: _D_<var>_<level>,
%                           numeric levels as integers or the shortest
%                           round-trip decimal (<= 17 digits).
%     fe               []   one or two fixed-effect columns (absorbed).
%     se         'robust'   'classical' | 'robust' (HC1) | 'cluster' |
%                           'spatial' (Conley).
%     clusterid        []   cluster column; used only with se='cluster'.
%     alpha          0.05   level of the Wald CIs and AR sets, in (0,1).
%     weights          []   weight column (non-negative; zero-weight rows
%                           are dropped).
%     wtype     'aweight'   'aweight' | 'pweight' (pweight + classical is
%                           upgraded to robust).
%     op_tau         0.10   Olea-Pflueger bias tolerance, in (0,1).
%     op_alpha       0.05   Olea-Pflueger test size, in (0,1).
%     ar_ci          true   Anderson-Rubin confidence sets.
%     ar_grid_halfwidth 200 K >= 2 search grid half-width (SE units).
%     ar_grid_points 4001   K >= 2 search grid points (integer >= 3).
%     lat, lon         []   coordinate columns (se='spatial' only).
%     spatial_cutoff   []   distance bandwidth (km under haversine).
%     spatial_kernel   []   'uniform' (default with se='spatial') |
%                           'bartlett'.
%     spatial_metric   []   'haversine' (default with se='spatial';
%                           degrees) | 'euclidean'.
%                           Spatial arguments with another se are an error.
%     lasso         false   post-double-selection LASSO (BCH 2014).
%     lasso_penalty 'plugin' 'plugin' only: 'cv' raises an error in this
%                           port (sklearn's folds cannot be reproduced).
%     lasso_keep       []   rich-expansion names always kept (globs OK).
%     lasso_seed 20260317   CV seed (unused by the plug-in penalty).
%     lasso_cv_folds    5   CV folds (validated; CV unavailable here).
%     expand_x {'square','cube'}  ops on continuous X: 'square', 'cube',
%                           'interact_with_X'.
%     expand_dummies {'interact_with_X'}  ops on dummy indicators:
%                           'interact_with_X', 'interact_with_dummies'.
%     expand_override  []   per-variable ops: containers.Map (name ->
%                           cellstr), an Nx2 cell {name, {ops}} or a
%                           struct (names such as '_D_v_2' need the Map
%                           or cell form).
%     show_selection false  print the LASSO selection details.
%     lasso_c_slack   1.1   plug-in constant c (> 0).
%     lasso_alpha_p    []   plug-in alpha_p, in (0,1); [] = 0.1/log(n).
%     lasso_n_iter_max 15   plug-in iterations (positive integer).
%     lasso_tol      1e-6   plug-in convergence tolerance (> 0).
%     lasso_lambda_override []  fixed penalty level lambda (>= 0).
%     lasso_verbose false   print the plug-in iterations.
%     quiet         false   suppress notes and warnings (MATLAB only).
%     ascii            []   plain-ASCII symbols (theta, gamma, ...)
%                           instead of Greek letters; [] = on in
%                           'matlab -batch' on Windows (whose console
%                           code page cannot show them), off otherwise.
%   Booleans are true/false or 0/1; text such as 'false' is an error.
%
%   RESULTS
%   -------
%   Returns a struct with the estimates, SEs, CIs, diagnostics (version,
%   n, n_clusters, pi1, se_pi1, F_TZ, F_MZ_T, *_op, K_eff_TZ, K_eff_MZ_T,
%   *_op_critval, AR sets, c_TM_hat, kappa_star, kappa_star_CI,
%   kappa_star_CI_informative, rho_TY_star, kappa_adj_table,
%   F_MZ_T_above_crit (true iff the effective F_{M,Z|T,X} reaches the
%   OP critical value; [] when either is not finite), F_reduced_form,
%   F_first_stage_cond, RV/XRV, df_fe_absorbed,
%   fe_nested_in_cluster, dropped_collinear, lasso_info, rich_col_names,
%   eps_*_hat, ...). Python None is [] here.
%
%   CONVENTIONS
%   -----------
%   - Wald CIs use normal critical values z_{1-alpha/2}; AR sets use chi2_K.
%   - Absorbed FE count as the rank of the FE dummies (G - 1 one-way;
%     G1 + G2 - c - 1 two-way, c connected components). With cluster SEs, FE
%     nested in the clusters are not counted in the CR1 correction
%     (reghdfe/fixest; with one of two dimensions nested, the other counts
%     G - 1); HC1 and classical SEs count all absorbed FE.
%   - Olea-Pflueger critical values: simplified test at the estimated K_eff
%     (Montiel Olea and Pflueger 2013, eq. 4), from a toolbox-free
%     non-central chi-square quantile; K_eff = 1 for one instrument
%     (23.109 at tau = 10%, alpha = 5%).
%   - RV/XRV with K >= 2 instruments use the first instrument.
%   - LASSO alpha_p = 0.1/log(n) in every regime (lassopack uses 0.1/log(G)
%     under clustering).
%
%   PORT NOTES (differences from the Python reference)
%   ---------------------------------------------------
%   - lasso_penalty='cv' is not available, and lasso_selection_stability
%     exists in the Python package only (their fold / subsample draws
%     cannot reproduce numpy's).
%   - Bootstrap draws differ from numpy's (the harness replays Python's).
%   - Integral numeric dummy levels are labeled '_D_v_2' (pandas float
%     columns give '_D_v_2.0'), as in the R package.
%   - Glob matching is case-sensitive (Python's fnmatch is case-insensitive
%     on Windows only).
%
%   REQUIREMENTS: MATLAB R2016b or later, no toolboxes (validated on R2026b).
%   GNU Octave: not tested; version 2.0.0 has not been run under Octave.
%
%   Version 2.0.0. Authors: Andreas Ferrara (maintainer, a.ferrara@pitt.edu),
%   Robert Gold, Stephan Heblich, Rodrigo Pinto.
%   Copyright (c) 2026 Andreas Ferrara, Robert Gold, Stephan Heblich, Rodrigo Pinto
%   License: MIT (see the LICENSE file).
%
%   Citation: Ferrara, Andreas, Robert Gold, Stephan Heblich, and Rodrigo
%   Pinto (2026), "Identifying Causal Mediation Effects with a Single
%   Instrument," NBER Working Paper No. [number], doi:[DOI].
%
%   See also IVMEDIATE_PRINT, IVMEDIATE_BOOTSTRAP, IVMEDIATE_PRINT_BOOTSTRAP,
%   IVMEDIATE_READ_CSV.

% Testing hook, not part of the interface: ivmediate('__internal__') returns
% the internal helpers of private/ivmediate_core.m (a struct of function
% handles) to the validation harness and the regression tests, since MATLAB
% does not allow private/ on the path.
if nargin == 1 && ischar(data) && strcmp(data, '__internal__')
    results = ivmediate_core('handles');
    return;
end
results = ivmediate_core('ivmediate', data, Y, T, M, Z, varargin{:});
end
