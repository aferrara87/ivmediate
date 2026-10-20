function ivmediate_print_bootstrap(boot, varargin)
%IVMEDIATE_PRINT_BOOTSTRAP  Print bootstrap mediation results (version 2.0.0).
%
%   ivmediate_print_bootstrap(boot)
%   ivmediate_print_bootstrap(boot, Name, Value, ...)
%
%   boot: the output of ivmediate_bootstrap. Options (names
%   case-insensitive; defaults in brackets): title ('Bootstrap Mediation'),
%   stars (true), star_cutoffs ([0.01 0.05 0.10]), ascii ([]: as in
%   ivmediate).
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
%   See also IVMEDIATE_BOOTSTRAP, IVMEDIATE_PRINT, IVMEDIATE.

ivmediate_core('ivmediate_print_bootstrap', boot, varargin{:});
end
