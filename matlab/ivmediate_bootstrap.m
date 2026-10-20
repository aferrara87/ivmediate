function boot = ivmediate_bootstrap(data, Y, T, M, Z, varargin)
%IVMEDIATE_BOOTSTRAP  Bootstrap inference for IV mediation (version 2.0.0).
%
%   boot = ivmediate_bootstrap(data, Y, T, M, Z)
%   boot = ivmediate_bootstrap(data, Y, T, M, Z, Name, Value, ...)
%
%   data, Y, T, M and Z as in ivmediate. Option names are case-insensitive;
%   defaults in brackets. Model options as in ivmediate (X, dummies, fe, se,
%   clusterid, alpha, weights, wtype, op_tau, op_alpha, lat, lon,
%   spatial_cutoff, spatial_kernel, spatial_metric; no LASSO), plus
%     method           []   'pairs' | 'block'; [] = 'block' when
%                           se='cluster', else 'pairs'.
%     B               999   replicates (integer >= 50).
%     seed       20260317   non-negative integer.
%     ci_type     'basic'   'basic' | 'percentile' | 'studentized'.
%     return_boot_dist false  add the replicate vectors (boot_dist).
%     bootstrap_diagnostics false  also bootstrap RV and XRV.
%     progress       true   progress lines.
%     verbose       false   report failed replicates.
%     quiet         false   suppress notes.
%     ascii            []   as in ivmediate.
%     boot_indices     []   testing hook: B x n_est (pairs) or B x G
%                           (block) 0-based draws replayed instead of
%                           the generator.
%     return_indices false  testing hook: out.indices (0-based draws)
%                           and out.rows (estimation-sample rows, 0-based).
%   The estimation sample is resampled (listwise deletion, zero weights
%   and, for block, missing cluster ids removed first). Draws come from a
%   portable splitmix64 counter generator keyed by a full-avalanche hash
%   of (seed, replicate); MATLAB's global random state is never used.
%   Basic/percentile replicates without diagnostics compute point
%   estimates only.
%
%   ivmediate_print_bootstrap(boot) prints the bootstrap results, and
%   ivmediate_print(results, 'boot', boot) puts the bootstrap SEs and CIs
%   in the coefficient table of the ivmediate results.
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
%   See also IVMEDIATE, IVMEDIATE_PRINT, IVMEDIATE_PRINT_BOOTSTRAP.

boot = ivmediate_core('ivmediate_bootstrap', data, Y, T, M, Z, varargin{:});
end
