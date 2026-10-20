% ivmediate: IV mediation analysis with a single instrument
% Version 2.0.0 17-Oct-2026
%
% Functions
%   ivmediate                 - Estimate the IV mediation model (results struct).
%   ivmediate_print           - Print the results of ivmediate.
%   ivmediate_bootstrap       - Bootstrap inference (pairs or cluster block).
%   ivmediate_print_bootstrap - Print the results of ivmediate_bootstrap.
%   ivmediate_read_csv        - Read a CSV file (NA and similar tokens as missing).
%
% Type "help ivmediate" for the options and the results. Example:
% examples/becker_woessmann/bw.m in the repository; smoke test:
% tests/test_bw.m (the paper's Table 6, column 1).
%
% Requires MATLAB R2016b or later, no toolboxes (validated on R2026b).
% GNU Octave: not tested.
%
% Authors: Andreas Ferrara (maintainer, a.ferrara@pitt.edu), Robert Gold,
% Stephan Heblich, Rodrigo Pinto.
% Copyright (c) 2026 Andreas Ferrara, Robert Gold, Stephan Heblich, Rodrigo Pinto
% License: MIT (see the LICENSE file).
%
% Citation: Ferrara, Andreas, Robert Gold, Stephan Heblich, and Rodrigo Pinto
% (2026), "Identifying Causal Mediation Effects with a Single Instrument,"
% NBER Working Paper No. [number], doi:[DOI].
