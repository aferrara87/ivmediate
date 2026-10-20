function ivmediate_print(results, varargin)
%IVMEDIATE_PRINT  Print IV mediation results (version 2.0.0).
%
%   ivmediate_print(results)
%   ivmediate_print(results, Name, Value, ...)
%
%   results: the output of ivmediate. Options (names case-insensitive;
%   defaults in brackets): title ('IV Mediation Analysis'), boot ([]:
%   bootstrap SEs/CIs in the table, from the output of ivmediate_bootstrap),
%   stars (true), star_cutoffs ([0.01 0.05 0.10]), ascii ([]: as in
%   ivmediate),
%     kappa_table   false   also print the kappa-adjusted estimates
%                           (theta and its CI with kappa_0 at 0, 10,
%                           25, 50, 75, 90 and 100% of kappa*); when
%                           F_MZ_T_above_crit is false, a note below
%                           the table says that c_TM is then estimated
%                           imprecisely.
%     r2_sensitivity false  also print the R^2-based sensitivity
%                           section (reduced-form F = t^2, RV, XRV and
%                           their components).
%   The kappa section always prints c_TM, kappa*, kappa*_CI and
%   rho*_TY (and the bootstrap CIs of kappa* and rho*_TY with 'boot');
%   every value is in the results whatever is printed.
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
%   See also IVMEDIATE, IVMEDIATE_BOOTSTRAP, IVMEDIATE_PRINT_BOOTSTRAP.

ivmediate_core('ivmediate_print', results, varargin{:});
end
