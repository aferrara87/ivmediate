function data = ivmediate_read_csv(path)
%IVMEDIATE_READ_CSV  Read a CSV file for ivmediate (version 2.0.0).
%
%   data = ivmediate_read_csv(path)
%
%   Reads a CSV file into a data set for ivmediate and ivmediate_bootstrap
%   (which also take a table or a struct of columns).
%
%   MATLAB: readtable with UTF-8 encoding; every column is read as text
%   and then typed: numeric when every non-missing field parses as a
%   number (Inf/-Inf included), text otherwise. Empty fields and the
%   tokens NA, N/A, NULL, NaN, nan and #N/A are missing in every column.
%   Quoted fields (RFC 4180) are supported. Octave (or a MATLAB whose
%   readtable lacks these options): a built-in parser with the same
%   typing and missing-value rules (UTF-8 BOM stripped, CRLF handled,
%   blank lines skipped; newlines inside quoted fields not supported).
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
%   See also IVMEDIATE, IVMEDIATE_BOOTSTRAP.

data = ivmediate_core('ivmediate_read_csv', path);
end
