function varargout = ivmediate_core(action, varargin)
%IVMEDIATE_CORE  Implementation of the ivmediate package (private).
%
%   MATLAB lets only the files of the parent folder call the functions in
%   private/, so this file is reached through the public functions there,
%   whose help text is the user documentation:
%       ivmediate(...)                  ivmediate_core('ivmediate', ...)
%       ivmediate_print(...)            ivmediate_core('ivmediate_print', ...)
%       ivmediate_bootstrap(...)        ivmediate_core('ivmediate_bootstrap', ...)
%       ivmediate_print_bootstrap(...)  ivmediate_core('ivmediate_print_bootstrap', ...)
%       ivmediate_read_csv(...)         ivmediate_core('ivmediate_read_csv', ...)
%
%   ivmediate_core('handles') returns a struct of function handles to the
%   internal helpers, and the version, for the validation harness and the
%   regression tests. They reach it through the testing hook
%   ivmediate('__internal__'), which is not part of the interface; a
%   function handle stays callable from any folder.
%   ivmediate_core('version') returns the version ('2.0.0').
%
%   The Python package (python/src/ivmediate/_core.py) is the reference
%   implementation; validation/run_matlab.m runs the cross-language harness.
%   Single file, base MATLAB only (no toolboxes). The Octave code paths are
%   kept (byte strings in the text helpers, the built-in CSV parser), but
%   version 2.0.0 has not been run under GNU Octave.
%
%   Version 2.0.0. Copyright (c) 2026 Andreas Ferrara, Robert Gold, Stephan
%   Heblich, Rodrigo Pinto. MIT License (see LICENSE in the parent folder).

if nargin < 1
    error('ivmediate_core: an action is required.');
end
switch action
    case 'ivmediate'
        varargout{1} = ivm_analysis(varargin{:});
    case 'ivmediate_print'
        ivm_print_results(varargin{:});
    case 'ivmediate_bootstrap'
        varargout{1} = ivm_bootstrap(varargin{:});
    case 'ivmediate_print_bootstrap'
        ivm_print_bootstrap(varargin{:});
    case 'ivmediate_read_csv'
        varargout{1} = ivm_read_csv(varargin{:});
    case 'handles'
        varargout{1} = ivm_handles();
    case 'version'
        varargout{1} = ivm_version();
    otherwise
        error('ivmediate_core: unknown action ''%s''.', action);
end
end

function v = ivm_version()
v = '2.0.0';
end


% ============================================================================
%  Internal helpers for the validation harness and the regression tests
% ============================================================================

function H = ivm_handles()
% Returned by the testing hook ivmediate('__internal__').
H = struct( ...
    'version', ivm_version(), ...
    'op_table', @ivm_op_table, ...
    'op_critical_value', @ivm_op_critical_value, ...
    'op_effective_F', @ivm_op_effective_F, ...
    'op_keff_critval', @ivm_op_keff_critval, ...
    'ncx2inv', @ivm_ncx2inv, ...
    'ncx2cdf', @ivm_ncx2cdf, ...
    'glob_to_regex', @ivm_glob_to_regex, ...
    'resolve_glob', @ivm_glob_filter_names, ...
    'resolve_var_list', @ivm_resolve_var_list, ...
    'prepare_sample', @ivm_prepare_sample, ...
    'norminv', @ivm_norminv, ...
    'normsf', @ivm_normsf, ...
    'tinv', @ivm_tinv, ...
    'chi2inv', @ivm_chi2inv, ...
    'chi2sf', @ivm_chi2sf, ...
    'chi2cdf', @ivm_chi2cdf, ...
    'quantile_linear', @ivm_quantile_linear, ...
    'sqrtnan', @ivm_sqrtnan, ...
    'eqqr', @ivm_eqqr, ...
    'lasso_cd', @ivm_lasso_cd, ...
    'bch_lasso', @ivm_bch_lasso, ...
    'compute_psi', @ivm_compute_psi, ...
    'spatial_adjacency', @ivm_spatial_adjacency, ...
    'twosls_with_vcov', @ivm_twosls_with_vcov, ...
    'ar_confidence_set', @ivm_ar_confidence_set, ...
    'weighted_cov', @ivm_wcov, ...
    'expand_dummies', @ivm_expand_dummies, ...
    'cat_label', @ivm_cat_label, ...
    'demean_twoway', @ivm_demean_twoway, ...
    'demean_within_groups', @ivm_demean_within_groups, ...
    'drop_singletons', @ivm_drop_singletons, ...
    'fe_df', @ivm_fe_df, ...
    'greedy_independent', @ivm_greedy_independent, ...
    'read_csv_custom', @ivm_read_csv_custom, ...
    'df_from_user', @ivm_df_from_user, ...
    'df_get', @ivm_df_get, ...
    'df_set', @ivm_df_set, ...
    'df_rows', @ivm_df_rows, ...
    'uch', @ivm_uch, ...
    'ascii_state', @ivm_ascii_state, ...
    'default_ascii', @ivm_default_ascii, ...
    'asciify', @ivm_asciify, ...
    'prng_key', @ivm_prng_key, ...
    'prng_uniform_int', @ivm_prng_uniform_int, ...
    'u64_mul', @ivm_u64_mul, ...
    'u64_mix', @ivm_u64_mix, ...
    'u64_from_hex', @ivm_u64_from_hex);
end

% ============================================================================
%  Environment, text output and option validation
% ============================================================================

function tf = ivm_is_octave()
persistent cached;
if isempty(cached)
    cached = exist('OCTAVE_VERSION', 'builtin') ~= 0;
end
tf = cached;
end

function s = ivm_uch_raw(cp)
%IVM_UCH_RAW Native string for a Unicode code point (BMP), never ASCII-mapped.
% Used for generated variable names (e.g., ' x ' with U+00D7 in interaction
% names), which must not depend on the display mode. MATLAB chars are UTF-16
% code units -> char(cp); Octave strings are bytes -> UTF-8 encoding.
if ~ivm_is_octave()
    s = char(cp);
    return;
end
if cp < 128
    s = char(cp);
elseif cp < 2048
    s = char([bitor(192, floor(cp / 64)), bitor(128, mod(cp, 64))]);
else
    s = char([bitor(224, floor(cp / 4096)), ...
              bitor(128, mod(floor(cp / 64), 64)), ...
              bitor(128, mod(cp, 64))]);
end
end

function s = ivm_uch(cp)
%IVM_UCH Display string for a Unicode code point: the native character, or
% its ASCII spelling while ASCII output is on (option 'ascii').
if ivm_ascii_state()
    s = ivm_ascii_name(cp);
else
    s = ivm_uch_raw(cp);
end
end

function s = ivm_ascii_name(cp)
% ASCII spelling of the symbols the printers use.
switch cp
    case 952,  s = 'theta';
    case 947,  s = 'gamma';
    case 964,  s = 'tau';
    case 954,  s = 'kappa';
    case 961,  s = 'rho';
    case 949,  s = 'eps';
    case 963,  s = 'sigma';
    case 960,  s = 'pi';
    case 955,  s = 'lambda';
    case 945,  s = 'alpha';
    case 946,  s = 'beta';
    case 916,  s = 'Delta';
    case 770,  s = '_hat';     % combining circumflex (theta-hat -> theta_hat)
    case 265,  s = 'c_hat';    % c with circumflex
    case 8320, s = '_0';
    case 8321, s = '_1';
    case 178,  s = '^2';
    case 8594, s = '->';
    case 9472, s = '-';
    case 8211, s = '-';
    case 8212, s = '--';
    case 215,  s = 'x';
    case 8746, s = 'U';
    case 8805, s = '>=';
    case 8804, s = '<=';
    otherwise
        if cp < 128
            s = char(cp);
        else
            s = '?';
        end
end
end

function s = ivm_asciify(s)
% Replace every non-ASCII character of s by its ASCII spelling.
if ivm_is_octave()
    b = double(s);
    if all(b < 128)
        return;
    end
    out = {};
    i = 1;
    while i <= numel(b)
        if b(i) < 128
            out{end+1} = char(b(i)); %#ok<AGROW>
            i = i + 1;
        elseif b(i) >= 224 && i + 2 <= numel(b)
            cp = (b(i) - 224) * 4096 + (b(i+1) - 128) * 64 + (b(i+2) - 128);
            out{end+1} = ivm_ascii_name(cp); %#ok<AGROW>
            i = i + 3;
        elseif b(i) >= 192 && i + 1 <= numel(b)
            cp = (b(i) - 192) * 64 + (b(i+1) - 128);
            out{end+1} = ivm_ascii_name(cp); %#ok<AGROW>
            i = i + 2;
        else
            out{end+1} = '?'; %#ok<AGROW>
            i = i + 1;
        end
    end
    s = [out{:}];
    return;
end
d = double(s);
k = find(d > 127);
if isempty(k)
    return;
end
parts = cell(1, numel(s));
for i = 1:numel(s)
    parts{i} = s(i);
end
for i = k
    parts{i} = ivm_ascii_name(d(i));
end
s = [parts{:}];
end

function tf = ivm_ascii_state(newval)
% Process-wide ASCII-output flag, set by the public functions for the
% duration of a call (option 'ascii').
persistent st;
if isempty(st)
    st = false;
end
if nargin > 0
    st = logical(newval);
end
tf = st;
end

function tf = ivm_resolve_ascii(v)
% 'ascii' option: [] means automatic (on in matlab -batch on Windows, whose
% console code page cannot show the Greek letters; off elsewhere).
if isempty(v) && ~ischar(v)
    tf = ivm_default_ascii();
else
    tf = ivm_check_bool('ascii', v);
end
end

function tf = ivm_default_ascii()
tf = false;
if ivm_is_octave() || ~ispc
    return;
end
try
    if exist('batchStartupOptionUsed', 'builtin') || exist('batchStartupOptionUsed', 'file')
        tf = logical(batchStartupOptionUsed);
    end
catch
    tf = false;
end
end

function ivm_out(fmt, varargin)
% All printed output goes through here (ASCII mode applies to every line).
s = sprintf(fmt, varargin{:});
if ivm_ascii_state()
    s = ivm_asciify(s);
end
fprintf(1, '%s', s);
end

function ivm_note(quiet, s)
% Print a note line unless quiet. s is a fully formatted string (no newline).
if ~quiet
    ivm_out('%s\n', s);
end
end

function n = ivm_strwidth(s)
%IVM_STRWIDTH Number of Unicode code points in s (display cells for padding),
% matching Python's str-length semantics used by format-spec padding.
if isempty(s)
    n = 0;
    return;
end
if ivm_is_octave()
    b = uint8(s);
    n = sum(bitand(b, uint8(192)) ~= uint8(128));   % count non-continuation bytes
else
    n = numel(s);
end
end

function out = ivm_pad_right(s, w)
% Left-align s in a field of width w (Python '{:<w}').
d = w - ivm_strwidth(s);
if d > 0
    out = [s, repmat(' ', 1, d)];
else
    out = s;
end
end

function out = ivm_pad_left(s, w)
% Right-align s in a field of width w (Python '{:>w}').
d = w - ivm_strwidth(s);
if d > 0
    out = [repmat(' ', 1, d), s];
else
    out = s;
end
end

function n = ivm_dwidth(s)
%IVM_DWIDTH Display width of s: code points without the combining
% diacritics U+0300..U+036F (the hat of theta-hat takes no column).
if isempty(s)
    n = 0;
    return;
end
if ivm_is_octave()
    b = double(uint8(s));
    n = sum(bitand(b, 192) ~= 128);
    for i = find((b == 204 | b == 205) & [true(1, numel(b) - 1), false])
        cp = (b(i) - 192) * 64 + (b(i + 1) - 128);
        n = n - (cp >= 768 && cp <= 879);
    end
else
    d = double(s);
    n = sum(~(d >= 768 & d <= 879));
end
end

function out = ivm_pad_right_dw(s, w)
% Left-align s in a field of w display columns (ivm_dwidth).
d = w - ivm_dwidth(s);
if d > 0
    out = [s, repmat(' ', 1, d)];
else
    out = s;
end
end

function w = ivm_label_width(labels, w0)
% Column width for a block of labels: Python's fixed width w0 in Unicode
% mode (output identical to the reference); in ASCII mode the longest ASCII
% label when it is wider, so every row still lines up with its header.
w = w0;
if ivm_ascii_state()
    for i = 1:numel(labels)
        w = max(w, ivm_strwidth(labels{i}));
    end
end
end

function s = ivm_thousands(n)
% Integer with thousands separators (Python '{:,}').
s = sprintf('%d', round(n));
neg = '';
if strncmp(s, '-', 1)
    neg = '-';
    s = s(2:end);
end
k = numel(s);
parts = {};
while k > 3
    parts = [{s(k-2:k)}, parts]; %#ok<AGROW>
    k = k - 3;
end
parts = [{s(1:k)}, parts];
s = [neg, strjoin(parts, ',')];
end

function c = ivm_norm_text_scalar(v, name)
% Normalize a scalar text value (char row / string scalar) to char.
if nargin < 2
    name = 'value';
end
if isa(v, 'string')
    if numel(v) ~= 1
        error('%s must be a single name; got a string array of %d elements.', name, numel(v));
    end
    c = char(v);
elseif ischar(v) && size(v, 1) <= 1
    c = v;
else
    error('%s must be text (a char vector or string scalar); got %s.', name, ivm_repr(v));
end
end

function out = ivm_norm_text_list(v, name)
% Normalize text-ish input to a cellstr row: char row -> {char},
% string array -> cellstr, cellstr -> cellstr. [] -> {}.
if nargin < 2
    name = 'value';
end
if isempty(v) && ~ischar(v)
    out = {};
elseif ischar(v)
    out = {v};
elseif isa(v, 'string')
    out = cellstr(v(:).');
elseif iscell(v)
    out = cell(1, numel(v));
    for i = 1:numel(v)
        out{i} = ivm_norm_text_scalar(v{i}, name);
    end
else
    error('%s must be text, a cell array of text, or a string array; got %s.', name, ivm_repr(v));
end
end

function s = ivm_repr(v)
% Short printable form of a value for error messages (Python repr-like).
if ischar(v)
    s = ['''', v, ''''];
elseif isa(v, 'string') && isscalar(v)
    s = ['"', char(v), '"'];
elseif islogical(v) && isscalar(v)
    if v
        s = 'true';
    else
        s = 'false';
    end
elseif isnumeric(v) && isscalar(v)
    s = sprintf('%.15g', v);
elseif isempty(v)
    s = '[]';
elseif isnumeric(v) || islogical(v)
    s = mat2str(v);
else
    s = sprintf('<%s>', class(v));
end
end

function s = ivm_pystr(v)
% Python str() of a float, approximately: integral floats keep '.0'.
if v == floor(v) && isfinite(v) && abs(v) < 1e16
    s = sprintf('%.1f', v);
else
    s = sprintf('%.15g', v);
end
end

function s = ivm_level_label(alpha)
% Confidence-level label used everywhere: 0.05 -> '95%', 0.025 -> '97.5%'.
s = sprintf('%g%%', round(100 * (1 - alpha) * 1e10) / 1e10);
end

function v = ivm_check_bool(name, v)
% Booleans must be logical scalars (true/false) or the numbers 0/1; text
% such as 'false' is an error.
if (islogical(v) || isnumeric(v)) && isscalar(v) && isreal(v) && (v == 0 || v == 1)
    v = logical(v);
    return;
end
error('%s must be true or false (a logical, or 0/1); got %s.', name, ivm_repr(v));
end

function v = ivm_check_choice(name, v, choices)
% Enumerated option: one of `choices` (exact, case-sensitive as in Python).
if isa(v, 'string') && isscalar(v)
    v = char(v);
end
if ~(ischar(v) && any(strcmp(choices, v)))
    q = cellfun(@(c) ['''', c, ''''], choices, 'UniformOutput', false);
    error('%s must be one of %s; got %s.', name, strjoin(q, ', '), ivm_repr(v));
end
end

function ivm_check_open_unit(name, v)
% alpha, op_tau, op_alpha and similar must lie strictly between 0 and 1.
if ~(isnumeric(v) && isscalar(v) && isreal(v))
    error('%s must be a number in (0, 1); got %s.', name, ivm_repr(v));
end
if ~(isfinite(v) && v > 0 && v < 1)
    error('%s must lie strictly between 0 and 1; got %s.', name, ivm_repr(v));
end
end

function v = ivm_check_posnum(name, v, allow_zero)
% A finite positive (or non-negative) number.
if nargin < 3
    allow_zero = false;
end
ok = isnumeric(v) && isscalar(v) && isreal(v) && isfinite(v);
if ok
    if allow_zero
        ok = v >= 0;
    else
        ok = v > 0;
    end
end
if ~ok
    if allow_zero
        error('%s must be a number >= 0; got %s.', name, ivm_repr(v));
    end
    error('%s must be a positive number; got %s.', name, ivm_repr(v));
end
v = double(v);
end

function v = ivm_check_int(name, v, lo)
% An integer-valued number >= lo.
if ~(isnumeric(v) && isscalar(v) && isreal(v) && isfinite(v) && v == round(v) && v >= lo)
    error('%s must be an integer >= %d; got %s.', name, lo, ivm_repr(v));
end
v = double(v);
end

function opts = ivm_parse_opts(defaults, args, fname)
% Name/value pairs; names are case-insensitive (MATLAB convention). Unknown
% names are an error that lists the valid ones.
opts = defaults;
if mod(numel(args), 2) ~= 0
    error('%s: options must be name/value pairs.', fname);
end
valid = fieldnames(defaults);
lv = lower(valid);
for i = 1:2:numel(args)
    key = args{i};
    if isa(key, 'string') && isscalar(key)
        key = char(key);
    end
    if ~ischar(key) || isempty(key)
        error('%s: option names must be text; got %s.', fname, ivm_repr(key));
    end
    hit = find(strcmp(lv, lower(key)), 1);
    if isempty(hit)
        error('%s: unknown option ''%s''. Valid options: %s.', fname, key, ...
              strjoin(valid(:).', ', '));
    end
    opts.(valid{hit}) = args{i + 1};
end
end

% ============================================================================
%  Internal dataframe: names (cellstr) + cols (cell of numeric col / cellstr)
%  Generated column names ('_D_x_a', 'x1^2', 'x1 <times> x2') are not valid
%  struct fields, so columns live in a parallel cell array.
% ============================================================================

function df = ivm_df_new(names, cols)
df = struct();
df.ivm_df_tag = true;
df.names = names(:).';
df.cols = cols(:).';
end

function tf = ivm_is_df(x)
tf = isstruct(x) && isscalar(x) && isfield(x, 'ivm_df_tag');
end

function names = ivm_data_names(data)
% Column names of a table, a struct of columns or an internal dataframe
% (no column is converted).
if ivm_is_df(data)
    names = data.names;
elseif isa(data, 'table')
    names = data.Properties.VariableNames;
elseif isstruct(data) && isscalar(data)
    names = fieldnames(data).';
else
    error(['data must be a table, a scalar struct of equal-length column ' ...
           'vectors, or the output of ivmediate_read_csv.']);
end
names = names(:).';
end

function df = ivm_df_from_user(data, cols)
%IVM_DF_FROM_USER Internal dataframe with the columns `cols` of data (all
% columns when omitted). Only these columns are converted, so unused columns
% of any type (datetime, duration, ...) are ignored.
names = ivm_data_names(data);
if nargin < 2
    cols = names;
end
c = cell(1, numel(cols));
for i = 1:numel(cols)
    nm = cols{i};
    ix = find(strcmp(names, nm), 1);
    if isempty(ix)
        error('column(s) not found in data: %s.', nm);
    end
    if ivm_is_df(data)
        c{i} = data.cols{ix};
    else
        c{i} = ivm_normalize_column(data.(nm), nm);
    end
end
df = ivm_df_new(cols, c);
ivm_df_check_lengths(df);
end

function col = ivm_normalize_column(v, name)
% A user column as double (numeric, logical) or cellstr (text types); text
% missing values become ''.
if isnumeric(v) || islogical(v)
    if ~isreal(v)
        error('Column ''%s'' is complex; use real numbers.', name);
    end
    col = double(v(:));
elseif iscell(v)
    if iscellstr(v)
        col = v(:);
        return;
    end
    isnum = cellfun(@(e) (isnumeric(e) || islogical(e)) && numel(e) <= 1, v(:));
    if all(isnum)
        col = NaN(numel(v), 1);
        for i = 1:numel(v)
            if ~isempty(v{i})
                col(i) = double(v{i});
            end
        end
        return;
    end
    col = cell(numel(v), 1);
    for i = 1:numel(v)
        e = v{i};
        if isempty(e)
            col{i} = '';
        elseif ischar(e) && size(e, 1) == 1
            col{i} = e;
        elseif isa(e, 'string') && isscalar(e)
            if ismissing(e)
                col{i} = '';
            else
                col{i} = char(e);
            end
        else
            error(['Column ''%s'' mixes text with other types (element %d is %s); ' ...
                   'use only text, or only numbers.'], name, i, class(e));
        end
    end
elseif isa(v, 'string')
    col = cell(numel(v), 1);
    miss = ismissing(v);
    for i = 1:numel(v)
        if miss(i)
            col{i} = '';
        else
            col{i} = char(v(i));
        end
    end
elseif isa(v, 'categorical')
    col = cellstr(v(:));
    col(isundefined(v(:))) = {''};
elseif ischar(v)
    col = cellstr(v);
else
    error('Column ''%s'' has unsupported type %s.', name, class(v));
end
end

function ivm_df_check_lengths(df)
if isempty(df.cols)
    return;
end
n = numel(df.cols{1});
for i = 2:numel(df.cols)
    if numel(df.cols{i}) ~= n
        error('Column ''%s'' has length %d; expected %d.', ...
              df.names{i}, numel(df.cols{i}), n);
    end
end
end

function n = ivm_df_nrows(df)
if isempty(df.cols)
    n = 0;
else
    n = numel(df.cols{1});
end
end

function tf = ivm_df_has(df, name)
tf = any(strcmp(df.names, name));
end

function col = ivm_df_get(df, name)
ix = find(strcmp(df.names, name), 1);
if isempty(ix)
    error('column ''%s'' not found in data.', name);
end
col = df.cols{ix};
end

function df = ivm_df_set(df, name, col)
ix = find(strcmp(df.names, name), 1);
if isempty(ix)
    df.names{end+1} = name;
    df.cols{end+1} = col;
else
    df.cols{ix} = col;
end
end

function df2 = ivm_df_rows(df, sel)
% Subset rows by logical mask or index vector.
df2 = df;
for i = 1:numel(df.cols)
    df2.cols{i} = df.cols{i}(sel);
end
end

function df2 = ivm_df_select(df, names)
% Keep the columns `names` (in that order).
c = cell(1, numel(names));
for i = 1:numel(names)
    c{i} = ivm_df_get(df, names{i});
end
df2 = ivm_df_new(names, c);
end

function m = ivm_col_isna(col)
if isnumeric(col)
    m = isnan(col);
else
    m = cellfun('isempty', col);
end
end

function x = ivm_col_num(col, name)
if ~isnumeric(col)
    error('column ''%s'' must be numeric for this role.', name);
end
x = double(col(:));
end

function [codes, G] = ivm_factorize(col)
% Integer codes 1..G in order of first appearance (pandas factorize(sort=False)).
if isnumeric(col)
    [~, ~, codes] = unique(col(:), 'stable');
else
    [~, ~, codes] = unique(col(:), 'stable');
end
codes = codes(:);
G = max([codes; 0]);
end


% ============================================================================
%  CSV reader
% ============================================================================

function tok = ivm_na_tokens()
% Tokens read as missing in any column (with the empty field).
tok = {'NA', 'N/A', 'NULL', 'NaN', 'nan', '#N/A'};
end

function df = ivm_read_csv(path)
%IVM_READ_CSV Load a CSV file into the internal dataframe (pass it to
% ivmediate / ivmediate_bootstrap like a table).
% MATLAB: readtable (UTF-8; quoted fields), every column read as text and
% typed here: a column is numeric when every non-missing token parses as a
% number ('Inf'/'-Inf' included), otherwise text. Empty fields and the tokens
% NA, N/A, NULL, NaN, nan and #N/A are missing in every column. Octave (or a
% MATLAB without these readtable options): the built-in parser below, with
% the same missing-value rules.
if isa(path, 'string') && isscalar(path)
    path = char(path);
end
if ~ischar(path)
    error('ivmediate_read_csv: path must be text; got %s.', ivm_repr(path));
end
if exist(path, 'file') ~= 2
    error('Cannot open file: %s', path);
end
df = [];
if ~ivm_is_octave()
    df = ivm_read_csv_readtable(path);
end
if isempty(df)
    df = ivm_read_csv_custom(path);
end
end

function df = ivm_read_csv_readtable(path)
% readtable with explicit options (no detectImportOptions heuristics, which
% can take a first data row with text tokens for a header line): the header
% is line 1 (parsed here), data lines 2..end, every column read as text.
df = [];
try
    fid = fopen(path, 'r', 'n', 'UTF-8');
    if fid < 0
        return;
    end
    header_line = fgetl(fid);
    fclose(fid);
    if ~ischar(header_line)
        error('CSV file %s is empty.', path);
    end
    if ~isempty(header_line) && double(header_line(1)) == 65279
        header_line = header_line(2:end);              % UTF-8 BOM
    end
    header_line = strrep(header_line, sprintf('\r'), '');
    header = ivm_split_csv_line(header_line);
    C = numel(header);
    opts = delimitedTextImportOptions('NumVariables', C, 'Delimiter', ',', ...
        'DataLines', [2, Inf], 'VariableNamesLine', 0, 'Encoding', 'UTF-8');
    opts.VariableNames = arrayfun(@(k) sprintf('v%d', k), 1:C, 'UniformOutput', false);
    opts.VariableTypes = repmat({'char'}, 1, C);
    opts.EmptyLineRule = 'skip';
    opts.ExtraColumnsRule = 'error';
    opts = setvaropts(opts, 'WhitespaceRule', 'preserve');
    T = readtable(path, opts);
catch err
    if strcmp(err.message, sprintf('CSV file %s is empty.', path))
        rethrow(err);
    end
    df = [];
    return;
end
nrows = height(T);
cols = cell(1, C);
for c = 1:C
    v = T.(sprintf('v%d', c));
    if ~iscell(v)
        v = cellstr(v);
    end
    cols{c} = v(:);
end
cols = ivm_type_columns(cols, nrows);
df = ivm_df_new(header, cols);
end

function df = ivm_read_csv_custom(path)
% Built-in CSV parser (used under Octave). Strips a UTF-8 BOM, handles CRLF,
% skips blank lines and handles double-quoted fields (RFC 4180; a quote
% character inside an unquoted field is literal). Newlines inside quoted
% fields are not supported.
fid = fopen(path, 'r');
if fid < 0
    error('Cannot open file: %s', path);
end
raw = fread(fid, Inf, 'uint8=>char').';
fclose(fid);
if numel(raw) >= 3 && double(raw(1)) == 239 && double(raw(2)) == 187 && double(raw(3)) == 191
    raw = raw(4:end);   % strip UTF-8 BOM
end
if ~ivm_is_octave()
    % MATLAB chars are UTF-16 code units: decode the bytes as UTF-8
    raw = native2unicode(uint8(raw), 'UTF-8');
end
raw = strrep(raw, sprintf('\r\n'), sprintf('\n'));
raw = strrep(raw, sprintf('\r'), sprintf('\n'));
lines = strsplit(raw, sprintf('\n'), 'CollapseDelimiters', false);
lines = lines(~cellfun(@(s) isempty(strtrim(s)), lines));
if isempty(lines)
    error('CSV file %s is empty.', path);
end
header = ivm_split_csv_line(lines{1});
C = numel(header);
nrows = numel(lines) - 1;
tokens = cell(nrows, C);
for i = 1:nrows
    f = ivm_split_csv_line(lines{i + 1});
    if numel(f) ~= C
        error('CSV row %d has %d fields; expected %d.', i, numel(f), C);
    end
    tokens(i, :) = f;
end
colcells = cell(1, C);
for c = 1:C
    colcells{c} = tokens(:, c);
end
cols = ivm_type_columns(colcells, nrows);
df = ivm_df_new(header, cols);
end

function cols = ivm_type_columns(colcells, nrows)
%IVM_TYPE_COLUMNS Decide numeric vs text per column and convert. Empty
% tokens and the NA tokens are missing; a column is numeric iff every other
% token is a decimal number ([+-] digits [. digits] [e[+-]digits], surrounding
% spaces allowed) or [+-]Inf / [+-]Infinity (any case), parsed exactly
% (sscanf, i.e., strtod: correctly rounded). Hexadecimal and other forms that
% strtod would accept ('0x1A', '1d3') are text, as in pandas.
na = ivm_na_tokens();
allowed = false(1, 128);
allowed(double('0123456789+-.eE ')) = true;
inf_words = {'inf', '+inf', '-inf', 'infinity', '+infinity', '-infinity'};
C = numel(colcells);
cols = cell(1, C);
for c = 1:C
    coltok = colcells{c}(:);
    miss = cellfun('isempty', coltok) | ismember(coltok, na);
    cj = coltok;
    cj(miss) = {'0'};
    if isempty(cj)
        cols{c} = zeros(0, 1);
        continue;
    end
    M = char(cj);                        % rows padded with spaces
    D = double(M);
    other = any(D >= 128, 2);
    Dk = D(~other, :);
    other(~other) = any(~allowed(Dk + (Dk == 0)), 2);   % code 0 -> index 1 (not allowed)
    if any(other)
        isinf_tok = ismember(lower(strtrim(cj(other))), inf_words);
        if ~all(isinf_tok)
            coltok(miss) = {''};
            cols{c} = coltok;
            continue;
        end
    end
    M(:, end + 1) = ' ';                 % guarantee a separator per field
    joined = reshape(M.', 1, []);
    [vals, cnt, ~, nextidx] = sscanf(joined, '%f');
    if cnt == nrows && (nextidx > numel(joined) || all(joined(nextidx:end) == ' '))
        vals = vals(:);
        vals(miss) = NaN;
        cols{c} = vals;
    else
        coltok(miss) = {''};
        cols{c} = coltok;
    end
end
end

function fields = ivm_split_csv_line(line)
if ~any(line == '"')
    fields = strsplit(line, ',', 'CollapseDelimiters', false);
    return;
end
fields = {};
cur = '';
in_q = false;
at_start = true;
i = 1;
L = numel(line);
while i <= L
    ch = line(i);
    if in_q
        if ch == '"'
            if i < L && line(i + 1) == '"'
                cur = [cur, '"']; %#ok<AGROW>
                i = i + 1;
            else
                in_q = false;
            end
        else
            cur = [cur, ch]; %#ok<AGROW>
        end
    else
        if ch == '"' && at_start
            in_q = true;
        elseif ch == ','
            fields{end+1} = cur; %#ok<AGROW>
            cur = '';
            at_start = true;
            i = i + 1;
            continue;
        else
            cur = [cur, ch]; %#ok<AGROW>
        end
    end
    at_start = false;
    i = i + 1;
end
fields{end+1} = cur;
end

% ============================================================================
%  Distribution functions (base numerics only; scipy-equivalent formulas)
% ============================================================================

function x = ivm_norminv(p)
x = -sqrt(2) .* erfcinv(2 .* p);
end

function p = ivm_normsf(x)
p = 0.5 .* erfc(x ./ sqrt(2));
end

function x = ivm_chi2inv(p, k)
x = 2 .* ivm_gammaincinv_c(p, k ./ 2);
end

function p = ivm_chi2sf(x, k)
p = gammainc(x ./ 2, k ./ 2, 'upper');
end

function p = ivm_chi2cdf(x, k)
p = gammainc(x ./ 2, k ./ 2);
end

function t = ivm_tinv(p, v)
%IVM_TINV Student-t quantile via the regularized incomplete beta inverse.
if p == 0.5
    t = 0.0;
    return;
end
if p < 0.5
    w = ivm_betaincinv_c(2 * p, v / 2, 0.5);
    t = -sqrt(v * (1 - w) / w);
else
    w = ivm_betaincinv_c(2 * (1 - p), v / 2, 0.5);
    t = sqrt(v * (1 - w) / w);
end
end

function x = ivm_gammaincinv_c(p, a)
%IVM_GAMMAINCINV_C Inverse regularized lower incomplete gamma. A built-in
% inverse (when available) is only a seed, always polished by Newton steps on
% gammainc, so the result is accurate in every part of the argument space.
if p <= 0
    x = 0; return;
end
if p >= 1
    x = Inf; return;
end
x = [];
if exist('gammaincinv', 'file') == 2 || exist('gammaincinv', 'builtin') ~= 0
    x = gammaincinv(p, a);
    if ~(isfinite(x) && x > 0)
        x = [];
    end
end
if isempty(x)
    z = ivm_norminv(p);                       % Wilson-Hilferty start
    x = a * (1 - 1 / (9 * a) + z * sqrt(1 / (9 * a)))^3;
    x = max(x, 1e-300);
end
lo = 0; hi = Inf;
for it = 1:200
    f = gammainc(x, a) - p;
    if f > 0
        hi = min(hi, x);
    else
        lo = max(lo, x);
    end
    d = exp((a - 1) * log(x) - x - gammaln(a));
    if d <= 0 || ~isfinite(d)
        xn = ivm_bisect_step(lo, hi);
    else
        xn = x - f / d;
        if ~(xn > lo && (isinf(hi) || xn < hi))
            xn = ivm_bisect_step(lo, hi);
        end
    end
    if abs(xn - x) <= 4 * eps * max(abs(x), 1e-300)
        x = xn;
        break;
    end
    x = xn;
end
end

function w = ivm_betaincinv_c(p, a, b)
%IVM_BETAINCINV_C Inverse regularized incomplete beta (built-in seed when
% available, then Newton on betainc; Octave's betaincinv loses digits at
% asymmetric arguments, which would corrupt small-alpha t quantiles).
if p <= 0
    w = 0; return;
end
if p >= 1
    w = 1; return;
end
w = [];
if exist('betaincinv', 'file') == 2 || exist('betaincinv', 'builtin') ~= 0
    w = betaincinv(p, a, b);
    if ~(isfinite(w) && w > 0 && w < 1)
        w = [];
    end
end
if isempty(w)
    w = min(max(p, 1e-12), 1 - 1e-12);
end
lo = 0; hi = 1;
lnB = gammaln(a) + gammaln(b) - gammaln(a + b);
for it = 1:200
    f = betainc(w, a, b) - p;
    if f > 0
        hi = min(hi, w);
    else
        lo = max(lo, w);
    end
    d = exp((a - 1) * log(w) + (b - 1) * log(1 - w) - lnB);
    if d <= 0 || ~isfinite(d)
        wn = 0.5 * (lo + hi);
    else
        wn = w - f / d;
        if ~(wn > lo && wn < hi)
            wn = 0.5 * (lo + hi);
        end
    end
    if abs(wn - w) <= 4 * eps * max(abs(w), 1e-300)
        w = wn;
        break;
    end
    w = wn;
end
end

function x = ivm_bisect_step(lo, hi)
if isinf(hi)
    x = max(2 * max(lo, 1), 1);
else
    x = 0.5 * (lo + hi);
end
end

function [w, j] = ivm_poisson_weights(h)
% Poisson(h) probabilities on the support that matters (mode +- 40 SDs),
% built by recursion from the mode and normalized to sum 1 (no gammaln).
jm = floor(h);
span = ceil(40 * sqrt(h)) + 40;
j = (max(0, jm - span):(jm + span)).';
w = zeros(size(j));
im = find(j == jm, 1);
w(im) = 1;
for i = im + 1:numel(j)
    w(i) = w(i - 1) * h / j(i);
end
for i = im - 1:-1:1
    w(i) = w(i + 1) * j(i + 1) / h;
end
w = w / sum(w);
end

function p = ivm_ncx2cdf(x, k, lam)
%IVM_NCX2CDF Noncentral chi-square CDF with df k > 0 (any real) and
% noncentrality lam >= 0: a Poisson(lam/2) mixture of central chi-square CDFs,
% sum_j Pois(j) * P(chi2_{k + 2j} <= x).
if x <= 0
    p = 0;
    return;
end
if lam <= 0
    p = gammainc(x / 2, k / 2);
    return;
end
[w, j] = ivm_poisson_weights(lam / 2);
p = sum(w .* gammainc(x / 2, k / 2 + j));
p = min(max(p, 0), 1);
end

function f = ivm_ncx2pdf(x, k, lam)
% Density of the same mixture.
if x <= 0
    f = 0;
    return;
end
if lam <= 0
    a = k / 2;
    f = exp((a - 1) * log(x / 2) - x / 2 - gammaln(a)) / 2;
    return;
end
[w, j] = ivm_poisson_weights(lam / 2);
a = k / 2 + j;
f = sum(w .* exp((a - 1) .* log(x / 2) - x / 2 - gammaln(a)) / 2);
end

function x = ivm_ncx2inv(p, k, lam)
%IVM_NCX2INV Noncentral chi-square quantile (toolbox-free): Brent's method
% (fzero) on the Poisson-mixture CDF, then Newton steps on the density.
if ~(p > 0 && p < 1) || ~(k > 0) || ~(lam >= 0) || ~isfinite(k) || ~isfinite(lam)
    if p == 0
        x = 0;
    elseif p == 1
        x = Inf;
    else
        x = NaN;
    end
    return;
end
F = @(t) ivm_ncx2cdf(t, k, lam) - p;
mu = k + lam;
sd = sqrt(2 * (k + 2 * lam));
hi = mu + 10 * sd + 10;
while F(hi) < 0
    hi = 2 * hi;
end
lo = 0;
if F(lo) > 0
    x = 0;
    return;
end
try
    x = fzero(F, [lo, hi], optimset('TolX', 1e-16 * hi, 'Display', 'off'));
catch
    x = NaN;
    return;
end
for it = 1:3
    fx = F(x);
    d = ivm_ncx2pdf(x, k, lam);
    if ~(d > 0) || ~isfinite(fx)
        break;
    end
    xn = x - fx / d;
    if ~(xn > 0) || abs(xn - x) <= 2 * eps * x
        if xn > 0
            x = xn;
        end
        break;
    end
    x = xn;
end
end


% ============================================================================
%  numpy-compat numeric utilities
% ============================================================================

function q = ivm_quantile_linear(x, p)
%IVM_QUANTILE_LINEAR numpy default quantile ('linear', a.k.a. R type 7).
xs = sort(x(:));
n = numel(xs);
q = zeros(size(p));
for i = 1:numel(p)
    if n == 0
        q(i) = NaN;
        continue;
    end
    if n == 1
        q(i) = xs(1);
        continue;
    end
    h = (n - 1) * p(i);
    lo = floor(h);
    frac = h - lo;
    hi = min(lo + 1, n - 1);
    q(i) = xs(lo + 1) * (1 - frac) + xs(hi + 1) * frac;
end
end

function s = ivm_sqrtnan(v)
% Real sqrt that propagates NaN and gives NaN (not complex) for negative
% inputs, as np.sqrt: a NaN or negative variance gives a NaN SE.
s = NaN(size(v));
ok = v >= 0;                 % false for NaN
s(ok) = sqrt(v(ok));
end

function tf = ivm_all_same_weight(w)
% np.allclose(w, w[0]) with numpy defaults rtol=1e-5, atol=1e-8.
tf = all(abs(w - w(1)) <= 1e-8 + 1e-5 .* abs(w(1)));
end


% ============================================================================
%  Bootstrap random numbers (portable, no global RNG state)
%
%  A splitmix64 counter generator on exact 64-bit integer arithmetic (four
%  16-bit limbs in doubles, vectorized). The key of replicate b is a
%  full-avalanche hash of (seed, b): two splitmix64 mixing rounds,
%  key_b = mix64(mix64(seed + gamma) + b * gamma). Draw i of replicate b is
%  mix64(key_b + i * gamma) (the splitmix64 stream started at key_b), so two
%  replicates share draws only if their keys differ by a multiple of gamma
%  smaller than the number of draws: probability about 2n/2^64 per pair.
%  MATLAB's global random-number state is never used or changed.
% ============================================================================

function L = ivm_u64_from_hex(h)
% 64-bit constant (16 hex digits) as 1x4 limbs, least significant first.
L = [hex2dec(h(13:16)), hex2dec(h(9:12)), hex2dec(h(5:8)), hex2dec(h(1:4))];
end

function L = ivm_u64_from_double(v)
% Non-negative integers below 2^53 (column vector) as n x 4 limbs.
v = double(v(:));
L = zeros(numel(v), 4);
for i = 1:4
    L(:, i) = mod(v, 65536);
    v = floor(v / 65536);
end
end

function L = ivm_u64_norm(S)
% Carry propagation of limb sums (each < 2^52) modulo 2^64.
L = S;
for i = 1:3
    c = floor(L(:, i) / 65536);
    L(:, i) = L(:, i) - c * 65536;
    L(:, i + 1) = L(:, i + 1) + c;
end
L(:, 4) = mod(L(:, 4), 65536);
end

function L = ivm_u64_add(A, B)
L = ivm_u64_norm(A + B);
end

function L = ivm_u64_mul(A, B)
% (A * B) mod 2^64; A, B are n x 4 or 1 x 4 limb arrays (broadcast).
S = zeros(max(size(A, 1), size(B, 1)), 4);
for i = 0:3
    for j = 0:(3 - i)
        S(:, i + j + 1) = S(:, i + j + 1) + A(:, i + 1) .* B(:, j + 1);
    end
end
% each product < 2^32 and at most 4 of them per limb: exact in doubles; split
% before normalizing so the carries stay exact
L = zeros(size(S));
carry = zeros(size(S, 1), 1);
for i = 1:4
    t = S(:, i) + carry;
    carry = floor(t / 65536);
    L(:, i) = t - carry * 65536;
end
end

function L = ivm_u64_xor(A, B)
L = bitxor(A, B);
end

function L = ivm_u64_shr(A, s)
% Logical right shift by s bits (0 < s < 64).
q = floor(s / 16);
r = s - 16 * q;
B = zeros(size(A));
B(:, 1:4 - q) = A(:, 1 + q:4);
if r == 0
    L = B;
    return;
end
L = floor(B / 2^r);
L(:, 1:3) = L(:, 1:3) + mod(B(:, 2:4), 2^r) * 2^(16 - r);
end

function z = ivm_u64_mix(z)
% splitmix64 finalizer.
persistent C1 C2
if isempty(C1)
    C1 = ivm_u64_from_hex('BF58476D1CE4E5B9');
    C2 = ivm_u64_from_hex('94D049BB133111EB');
end
z = ivm_u64_mul(ivm_u64_xor(z, ivm_u64_shr(z, 30)), C1);
z = ivm_u64_mul(ivm_u64_xor(z, ivm_u64_shr(z, 27)), C2);
z = ivm_u64_xor(z, ivm_u64_shr(z, 31));
end

function g = ivm_u64_gamma()
g = ivm_u64_from_hex('9E3779B97F4A7C15');
end

function key = ivm_prng_key(seed, b)
% Full-avalanche key of replicate b (1-based) for a non-negative integer seed.
g = ivm_u64_gamma();
s = ivm_u64_mix(ivm_u64_add(ivm_u64_from_double(seed), g));
key = ivm_u64_mix(ivm_u64_add(s, ivm_u64_mul(ivm_u64_from_double(b), g)));
end

function idx = ivm_prng_uniform_int(seed, b, count, n)
%IVM_PRNG_UNIFORM_INT `count` draws uniform on 1..n for replicate b of `seed`
% (count x 1). Draw i is mix64(key_b + i*gamma); its top 53 bits give
% u in [0, 1) and the draw is 1 + floor(u * n).
key = ivm_prng_key(seed, b);
st = ivm_u64_add(repmat(key, count, 1), ivm_u64_mul(ivm_u64_from_double((1:count).'), ivm_u64_gamma()));
z = ivm_u64_mix(st);
hi = z(:, 4) * 65536 + z(:, 3);                      % top 32 bits
lo = z(:, 2) * 65536 + z(:, 1);                      % low 32 bits
u = (hi * 2^21 + floor(lo / 2^11)) / 2^53;           % 53-bit uniform in [0, 1)
idx = 1 + floor(u * n);
idx(idx > n) = n;
end

% ============================================================================
%  Hajek-form weighted moments
% ============================================================================

function m = ivm_wmean(x, w)
m = sum(w .* x) / sum(w);
end

function c = ivm_wcov(x, y, w)
wsum = sum(w);
mx = sum(w .* x) / wsum;
my = sum(w .* y) / wsum;
c = sum(w .* (x - mx) .* (y - my)) / wsum;
end

function v = ivm_wvar(x, w)
v = ivm_wcov(x, x, w);
end

function s = ivm_wstd(x, w)
s = sqrt(ivm_wvar(x, w));
end


% ============================================================================
%  Least squares (2026-10: QR of the equilibrated design, never X'X)
%
%  With A*diag(scale) = Q*R (columns scaled to unit norm), the 2SLS bread is
%  (Xh'Xh)^{-1} = T*T' with T = diag(scale)*R^{-1}, and the score of row i in
%  the rotated coordinates is Q(i,:)*e(i), so every sandwich is
%  V = T*meat_Q*T'. Ill-conditioned but full-rank designs (calendar-year
%  polynomials, variables in large units) stay accurate.
% ============================================================================

function [Q, R, scale] = ivm_eqqr(A)
%IVM_EQQR Economy QR of A with columns scaled to unit Euclidean norm; error
% if the scaled matrix is numerically singular (|R_jj| <= 1e-12*max|R_ii|).
norms = sqrt(sum(A .* A, 1));
scale = ones(1, size(A, 2));
nz = norms > 0;
scale(nz) = 1 ./ norms(nz);
[Q, R] = qr(A .* scale, 0);
d = abs(diag(R));
if ~isempty(d) && (any(~isfinite(d)) || min(d) <= 1e-12 * max(max(d), 1e-300))
    error(['Singular matrix: the design matrix is numerically rank-deficient ' ...
           '(collinear columns).']);
end
scale = scale(:);
end

function fit = ivm_lm_fit(X, Y)
% OLS (np.linalg.lstsq): economy QR for full-rank designs, pinv (minimum-norm
% solution) for rank-deficient ones.
[Q, R] = qr(X, 0);
d = abs(diag(R));
if isempty(d) || min(d) <= max(size(X)) * eps(max(d))
    beta = pinv(X) * Y;
else
    beta = R \ (Q' * Y);
end
fit = struct();
fit.coefficients = beta;
fit.residuals = Y - X * beta;
end

function fit = ivm_wls_fit(X, Y, w)
% Weighted OLS on raw data via QR of sqrt(w)*X; raw-scale residuals.
sw = sqrt(w);
[Q, R, scale] = ivm_eqqr(X .* sw);
beta = scale .* (R \ (Q' * (Y .* sw)));
fit = struct();
fit.coefficients = beta;
fit.residuals = Y - X * beta;
end

function r = ivm_residualize(A, B)
% Residuals of the columns of A on B (QR of the equilibrated B).
if isempty(B) || size(B, 2) == 0
    r = A;
    return;
end
[Q, ~, ~] = ivm_eqqr(B);
r = A - Q * (Q' * A);
end

function tf = ivm_is_pd(V)
% True if V is finite and positive definite (a usable variance block).
if isempty(V) || any(~isfinite(V(:)))
    tf = false;
    return;
end
tf = min(eig((V + V') / 2)) > 0;
end


% ============================================================================
%  Olea-Pflueger (2013) effective F and critical values
% ============================================================================

function [tab, taus] = ivm_op_table()
% Olea & Pflueger (2013, JBES) Table 1 (alpha = 0.05), K = 1..30, for
% reference only: since 2026-10 critical values are computed at K_eff.
taus = [0.05, 0.10, 0.20, 0.30];
tab = [ ...
    37.42, 23.11, 15.06, 12.05;
    32.32, 19.29, 12.17,  9.57;
    30.13, 17.67, 10.95,  8.53;
    28.85, 16.72, 10.23,  7.92;
    27.98, 16.08,  9.75,  7.51;
    27.35, 15.62,  9.40,  7.21;
    26.86, 15.26,  9.14,  6.98;
    26.47, 14.97,  8.92,  6.80;
    26.15, 14.73,  8.74,  6.65;
    25.87, 14.53,  8.59,  6.52;
    25.64, 14.36,  8.47,  6.41;
    25.44, 14.21,  8.36,  6.32;
    25.26, 14.08,  8.26,  6.24;
    25.10, 13.96,  8.17,  6.16;
    24.96, 13.86,  8.10,  6.10;
    24.83, 13.77,  8.03,  6.04;
    24.71, 13.68,  7.96,  5.99;
    24.60, 13.60,  7.91,  5.94;
    24.50, 13.53,  7.85,  5.89;
    24.41, 13.46,  7.80,  5.85;
    24.33, 13.40,  7.76,  5.81;
    24.25, 13.35,  7.72,  5.78;
    24.18, 13.29,  7.68,  5.74;
    24.11, 13.24,  7.64,  5.71;
    24.05, 13.20,  7.61,  5.68;
    23.98, 13.15,  7.57,  5.66;
    23.93, 13.11,  7.54,  5.63;
    23.87, 13.07,  7.51,  5.61;
    23.82, 13.04,  7.49,  5.58;
    23.77, 13.00,  7.46,  5.56];
if abs(tab(1, 2) - 23.11) >= 1e-9
    error('OP table mis-keyed: K=1, tau=0.10 should be 23.11.');
end
end

function cv = ivm_op_critical_value(K, tau, alpha)
% Simplified-test critical value at K_eff = K: ncx2inv(1 - alpha, K, K/tau)/K
% (reproduces OP Table 1 to its two decimals at alpha = 0.05).
ivm_check_open_unit('tau', tau);
ivm_check_open_unit('alpha', alpha);
x = 1 / tau;
cv = ivm_ncx2inv(1 - alpha, K, x * K) / K;
end

function F = ivm_op_effective_F(beta_Z, V_ZZ, Q_zz)
beta_Z = beta_Z(:);
F = (beta_Z' * Q_zz * beta_Z) / trace(Q_zz * V_ZZ);
end

function [keff, crit] = ivm_op_keff_critval(V, Q, tau, alpha)
%IVM_OP_KEFF_CRITVAL Montiel Olea-Pflueger (2013) simplified test at the
% estimated K_eff: W2 = Q^{1/2} V Q^{1/2}; K_eff = [tr W2]^2 (1 + 2x) /
% (tr(W2'W2) + 2x tr(W2) lambda_max(W2)), x = 1/tau; critical value
% ncx2inv(1 - alpha, K_eff, x K_eff) / K_eff. K = 1 gives K_eff = 1.
x = 1 / tau;
K = size(V, 1);
if K == 1
    keff = 1;
else
    if ~(ivm_is_pd(V) && all(isfinite(Q(:))))
        keff = NaN;
        crit = NaN;
        return;
    end
    [U, ev] = eig((Q + Q') / 2);
    ev = diag(ev);
    Qh = (U .* sqrt(max(ev, 0)).') * U';
    W2 = Qh * V * Qh;
    W2 = (W2 + W2') / 2;
    tr = trace(W2);
    lmax = max(eig(W2));
    den = trace(W2 * W2) + 2 * x * tr * lmax;
    if ~(isfinite(den) && den > 0 && tr > 0)
        keff = NaN;
        crit = NaN;
        return;
    end
    keff = tr * tr * (1 + 2 * x) / den;
end
crit = ivm_ncx2inv(1 - alpha, keff, x * keff) / keff;
end


% ============================================================================
%  Variance-covariance estimators for 2SLS (rotated coordinates)
% ============================================================================

function vcov = ivm_vcov_classical(resid, T, n, k, df_adjust)
sigma2 = sum(resid .^ 2) / (n - k - df_adjust);
vcov = sigma2 * (T * T');
end

function vcov = ivm_vcov_robust(resid, Q, T, n, k, df_adjust)
% HC1: T*[(Q.*e)'(Q.*e)]*T' * n/(n-k-df).
S = Q .* resid;
meat = (S' * S) * (n / (n - k - df_adjust));
vcov = T * meat * T';
end

function vcov = ivm_vcov_cluster(resid, Q, T, n, k, gidx, G, df_adjust)
% CR1 (Liang-Zeger): correction G/(G-1)*(n-1)/(n-k-df); df excludes FE
% nested in the clusters (the caller passes that count).
if G < 2
    error('cluster SEs need at least 2 clusters.');
end
S = Q .* resid;
sums = sparse(gidx, 1:n, 1, G, n) * S;
correction = (G / (G - 1)) * ((n - 1) / (n - k - df_adjust));
meat = full(sums' * sums) * correction;
vcov = T * meat * T';
end

function [vcov, n_bad] = ivm_vcov_spatial(resid, Q, T, n, k, A, df_adjust)
% Conley: T*[scores'*A*scores]*T' * (n-1)/max(n-k-df, 1). The kernel matrix
% need not be PSD; negative variances are counted by the caller.
scores = Q .* resid;
meat = full(scores' * (A * scores)) * ((n - 1) / max(n - k - df_adjust, 1));
vcov = T * meat * T';
n_bad = sum(diag(vcov) < 0);
end

function [A, avg_neighbors, nnz_count] = ivm_spatial_adjacency(lat, lon, cutoff, kernel, metric)
%IVM_SPATIAL_ADJACENCY Kernel-weighted sparse adjacency A(i,j) = K(d_ij) for
% all pairs with d_ij <= cutoff (self included, d = 0). Haversine distance in
% km (Earth radius 6371) for lat/lon in degrees, or raw Euclidean distance.
% Distances are computed in row blocks to bound the temporary dense matrix.
EARTH_RADIUS_KM = 6371.0;
n = numel(lat);
is_hav = strcmp(metric, 'haversine');
if is_hav
    latr = lat(:) * (pi / 180);
    lonr = lon(:) * (pi / 180);
    coslat = cos(latr);
else
    px = lat(:);
    py = lon(:);
end
block = max(1, floor(5e6 / max(n, 1)));
ii = cell(1, ceil(n / block));
jj = ii;
kk = ii;
nb = 0;
nnz_count = 0;
start = 1;
while start <= n
    idx = start:min(start + block - 1, n);
    if is_hav
        dlat = latr(idx) - latr';
        dlon = lonr(idx) - lonr';
        a = sin(dlat / 2).^2 + (coslat(idx) * coslat') .* sin(dlon / 2).^2;
        a(a < 0) = 0;
        a(a > 1) = 1;
        D = 2 * EARTH_RADIUS_KM * asin(sqrt(a));
    else
        dx = px(idx) - px';
        dy = py(idx) - py';
        D = sqrt(dx .* dx + dy .* dy);
    end
    hit = find(D <= cutoff);
    nb = nb + 1;
    if ~isempty(hit)
        [r_loc, c_glob] = ind2sub(size(D), hit);
        d = D(hit);
        if strcmp(kernel, 'uniform')
            w = ones(numel(d), 1);
        else
            w = max(1.0 - d(:) / cutoff, 0.0);
        end
        ii{nb} = reshape(idx(r_loc), [], 1);
        jj{nb} = c_glob(:);
        kk{nb} = w(:);
        nnz_count = nnz_count + numel(d);
    end
    start = start + block;
end
A = sparse(vertcat(ii{:}), vertcat(jj{:}), vertcat(kk{:}), n, n);
avg_neighbors = nnz_count / max(n, 1);   % within-cutoff pairs incl. self
end

function cov = ivm_cross_equation_cov(res1, res2, idx1, idx2, sekw, df_adjust)
%IVM_CROSS_EQUATION_COV Cov(beta1_idx1, beta2_idx2) from a stacked-GMM
% treatment of two 2SLS regressions sharing observations (delta-method SE of
% the indirect effect). Influence functions psi = (Q*T(idx,:)') .* resid.
n = numel(res1.resid);
k_eff = max(size(res1.T, 1), size(res2.T, 1));
psi1 = (res1.Q * res1.T(idx1, :)') .* res1.resid;
psi2 = (res2.Q * res2.T(idx2, :)') .* res2.resid;
switch sekw.se
    case 'none'
        cov = NaN;
    case 'classical'
        sigma_12 = (res1.resid' * res2.resid) / max(n - k_eff - df_adjust, 1);
        cross_block = res1.T * (res1.Q' * res2.Q) * res2.T';
        cov = sigma_12 * cross_block(idx1, idx2);
    case 'robust'
        cov = sum(psi1 .* psi2) * (n / (n - k_eff - df_adjust));
    case 'cluster'
        G = sekw.cluster_G;
        s1 = accumarray(sekw.cluster_gidx, psi1, [G, 1]);
        s2 = accumarray(sekw.cluster_gidx, psi2, [G, 1]);
        cov = sum(s1 .* s2) * ((G / (G - 1)) * ((n - 1) / (n - k_eff - df_adjust)));
    case 'spatial'
        cov = (psi1' * (sekw.spatial_A * psi2)) * ((n - 1) / max(n - k_eff - df_adjust, 1));
    otherwise
        error('Unknown se type: %s', sekw.se);
end
end

function res = ivm_twosls_with_vcov(Y, X, Z, sekw, df_adjust)
%IVM_TWOSLS_WITH_VCOV 2SLS (OLS when Z = []) via QR of the equilibrated
% designs: Xh = P_Z X from a QR of Z, beta from a QR of Xh. Returns coef, se,
% vcov, resid, sigma2, df and the factors Q (Xh*T = Q) and T
% ((Xh'Xh)^{-1} = T*T'). sekw.se = 'none' skips the variance (NaN).
n = numel(Y);
k = size(X, 2);
effective_df = n - k - df_adjust;
if isempty(Z)
    Xh = X;
else
    [Qz, ~, ~] = ivm_eqqr(Z);
    Xh = Qz * (Qz' * X);
end
[Q, R, sc] = ivm_eqqr(Xh);
beta = sc .* (R \ (Q' * Y));
resid = Y - X * beta;
Tm = sc .* (R \ eye(k));
n_bad = 0;
switch sekw.se
    case 'none'
        vcov = NaN(k, k);
    case 'classical'
        vcov = ivm_vcov_classical(resid, Tm, n, k, df_adjust);
    case 'robust'
        vcov = ivm_vcov_robust(resid, Q, Tm, n, k, df_adjust);
    case 'cluster'
        vcov = ivm_vcov_cluster(resid, Q, Tm, n, k, sekw.cluster_gidx, ...
                                sekw.cluster_G, df_adjust);
    case 'spatial'
        [vcov, n_bad] = ivm_vcov_spatial(resid, Q, Tm, n, k, sekw.spatial_A, df_adjust);
    otherwise
        error(['Unknown se type: %s. Use ''classical'', ''robust'', ' ...
               '''cluster'', or ''spatial''.'], sekw.se);
end
res = struct();
res.coef = beta;
res.se = ivm_sqrtnan(diag(vcov));
res.vcov = vcov;
res.resid = resid;
res.sigma2 = sum(resid .^ 2) / effective_df;
res.df = effective_df;
res.Q = Q;
res.T = Tm;
res.n_bad = n_bad;
end

% ============================================================================
%  Anderson-Rubin confidence sets (weak-instrument-robust inference)
% ============================================================================

function out = ivm_ar_confidence_set(y_dep, x_endog, W_incl, z_cols, sekw, ...
                                     df_adjust, alpha, center, scale, ...
                                     grid_halfwidth, grid_points)
%IVM_AR_CONFIDENCE_SET AR set for the coefficient on one endogenous
% regressor, by test inversion under the active SE regime (Python reference
% _ar_confidence_set). The AR regression coefficients are linear and their
% variance block quadratic in b0 for every regime, so three regressions give
% the statistic in closed form.
%   K = 1: exact Fieller inversion (checked at a fourth point).
%   K > 1: closed-form statistic on a search grid (the grid only brackets
%          sign changes), boundaries refined by Brent's method, intervals
%          narrower than the grid step found at interior local minima.
%   Spatial SEs: a variance that turns negative within center +-
%   grid_halfwidth*scale makes the set 'unavailable' (NaN bounds).
K = numel(z_cols);
cv = ivm_chi2inv(1.0 - alpha, K);

[bz0, Vz0] = ivm_ar_eval(0.0, y_dep, x_endog, W_incl, z_cols, sekw, df_adjust);
if ivm_is_pd(Vz0)
    stat0 = bz0' * (Vz0 \ bz0);
    pval0 = ivm_chi2sf(stat0, K);
else
    stat0 = NaN;
    pval0 = NaN;
end

out = struct();
out.shape = '';
out.lower = NaN; out.upper = NaN;
out.left_upper = NaN; out.right_lower = NaN;
out.segments = struct('lo', {}, 'hi', {}, 'open_left', {}, 'open_right', {});
out.has_segments = false;
out.cv = cv; out.df_chi2 = K; out.alpha = alpha;
out.stat_at_zero = stat0; out.pval_at_zero = pval0;
out.contains_zero = [];
out.method = '';

if isfinite(center)
    tc = center;
else
    tc = 0.0;
end
if ~isempty(scale) && isfinite(scale) && scale > 0
    sc = scale;
else
    sc = max(abs(tc), 1.0);
end
s_hw = grid_halfwidth * sc;
check_neg = strcmp(sekw.se, 'spatial');
use_grid = (K > 1);

if ~use_grid
    % -- Exact Fieller inversion (K = 1) --
    d = sc;
    [b_m, V_m] = ivm_ar_eval(tc - d, y_dep, x_endog, W_incl, z_cols, sekw, df_adjust);
    [b_c, V_c] = ivm_ar_eval(tc,     y_dep, x_endog, W_incl, z_cols, sekw, df_adjust);
    [b_p, V_p] = ivm_ar_eval(tc + d, y_dep, x_endog, W_incl, z_cols, sekw, df_adjust);
    bm = b_m(1); b0c = b_c(1); bp = b_p(1);
    Vm = V_m(1, 1); V0 = V_c(1, 1); Vp = V_p(1, 1);
    B = (bm - bp) / (2.0 * d);
    bv = (Vp - Vm) / (2.0 * d);
    cvv = (Vp - 2.0 * V0 + Vm) / (2.0 * d * d);
    qa = B * B - cv * cvv;
    qb = -2.0 * b0c * B - cv * bv;
    qc = b0c * b0c - cv * V0;
    [b_chk, V_chk] = ivm_ar_eval(tc + 2.0 * d, y_dep, x_endog, W_incl, z_cols, sekw, df_adjust);
    q_pred = qa * (2.0 * d)^2 + qb * (2.0 * d) + qc;
    q_act = b_chk(1)^2 - cv * V_chk(1, 1);
    q_scale = max([abs(qa) * (2 * d)^2, abs(qb) * (2 * d), abs(qc), abs(q_act), 1e-300]);
    if ~isfinite(q_pred) || abs(q_pred - q_act) > 1e-6 * q_scale
        use_grid = true;   % defensive fallback; should not trigger
    else
        out.method = 'fieller';
        s_c = [-s_hw, 0.0, s_hw];
        if cvv > 0 && abs(bv / (2.0 * cvv)) <= s_hw
            s_c(end+1) = -bv / (2.0 * cvv);
        end
        v_neg = false;
        if check_neg
            for s = s_c
                if V0 + bv * s + cvv * s * s < -1e-12 * (abs(V0) + abs(bv * s) + abs(cvv) * s * s)
                    v_neg = true;
                end
            end
        end
        tiny = 1e-14 * max([abs(qb) * d, abs(qc), 1e-300]);
        if v_neg
            out.shape = 'unavailable';
        elseif abs(qa) * d * d <= tiny
            if abs(qb) * d <= tiny
                if qc <= 0
                    out.shape = 'line';
                else
                    out.shape = 'empty';
                end
            elseif qb > 0
                out.shape = 'rays_left_only';
                out.left_upper = tc + (-qc / qb);
            else
                out.shape = 'rays_right_only';
                out.right_lower = tc + (-qc / qb);
            end
        else
            disc = qb * qb - 4.0 * qa * qc;
            if disc < 0.0
                if qa > 0
                    out.shape = 'empty';
                else
                    out.shape = 'line';
                end
            else
                sq = sqrt(disc);
                if qb ~= 0
                    if qb < 0
                        q_ = -0.5 * (qb - sq);
                    else
                        q_ = -0.5 * (qb + sq);
                    end
                else
                    q_ = 0.5 * sq;
                end
                r1 = q_ / qa;
                if q_ ~= 0
                    r2 = qc / q_;
                else
                    r2 = -qb / qa - r1;
                end
                s_lo = min(r1, r2);
                s_hi = max(r1, r2);
                if qa > 0
                    out.shape = 'interval';
                    out.lower = tc + s_lo;
                    out.upper = tc + s_hi;
                else
                    out.shape = 'rays';
                    out.left_upper = tc + s_lo;
                    out.right_lower = tc + s_hi;
                end
            end
        end
    end
end

if use_grid
    % -- Closed-form statistic on a grid, boundaries refined (K > 1) --
    out.method = 'grid';
    out.has_segments = true;
    grid = tc + linspace(-grid_halfwidth, grid_halfwidth, round(grid_points)) * sc;
    n_g = numel(grid);
    d = sc;
    [b_m, V_m] = ivm_ar_eval(tc - d, y_dep, x_endog, W_incl, z_cols, sekw, df_adjust);
    [b_c, V_c] = ivm_ar_eval(tc,     y_dep, x_endog, W_incl, z_cols, sekw, df_adjust);
    [b_p, V_p] = ivm_ar_eval(tc + d, y_dep, x_endog, W_incl, z_cols, sekw, df_adjust);
    cf = struct('tc', tc, 'b_c', b_c, 'Bv', (b_m - b_p) / (2.0 * d), 'V_c', V_c, ...
                'V1', (V_p - V_m) / (2.0 * d), 'V2', (V_p - 2.0 * V_c + V_m) / (2.0 * d * d));
    brute = @(b0) ivm_ar_stat(b0, y_dep, x_endog, W_incl, z_cols, sekw, df_adjust);
    chk_b = tc + 2.0 * d;
    chk_act = brute(chk_b);
    chk_pred = ivm_ar_closed(cf, chk_b);
    exact = all(isfinite([chk_act, chk_pred])) && abs(chk_pred - chk_act) <= 1e-6 * max(abs(chk_act), 1.0);
    stat_grid = zeros(1, n_g);
    lam_lo = zeros(1, n_g);
    lam_sc = zeros(1, n_g);
    for i = 1:n_g
        if exact
            [stat_grid(i), ev] = ivm_ar_closed(cf, grid(i));
        else
            [bz, Vz] = ivm_ar_eval(grid(i), y_dep, x_endog, W_incl, z_cols, sekw, df_adjust);
            ev = eig((Vz + Vz') / 2);
            if ivm_is_pd(Vz)
                stat_grid(i) = bz' * (Vz \ bz);
            else
                stat_grid(i) = NaN;
            end
        end
        if check_neg
            lam_lo(i) = min(ev);
            lam_sc(i) = max(abs(ev));
        end
    end
    if exact
        f_ = @(b) ivm_ar_closed(cf, b) - cv;
    else
        f_ = @(b) brute(b) - cv;
    end
    if check_neg
        unavailable = any(~isfinite(lam_lo)) || any(lam_lo < -1e-12 * lam_sc);
    else
        unavailable = false;
    end
    accept = stat_grid <= cv;          % NaN -> false
    if unavailable
        out.shape = 'unavailable';
        out.segments = struct('lo', {}, 'hi', {}, 'open_left', {}, 'open_right', {});
        out.has_segments = false;
    else
        segs = struct('lo', {}, 'hi', {}, 'open_left', {}, 'open_right', {});
        i = 1;
        while i <= n_g
            if accept(i)
                j = i;
                while j + 1 <= n_g && accept(j + 1)
                    j = j + 1;
                end
                lo = grid(i);
                hi = grid(j);
                if i > 1
                    r_ = ivm_ar_root(f_, grid(i - 1), grid(i), sc);
                    if ~isempty(r_), lo = r_; end
                end
                if j < n_g
                    r_ = ivm_ar_root(f_, grid(j), grid(j + 1), sc);
                    if ~isempty(r_), hi = r_; end
                end
                segs(end+1) = struct('lo', lo, 'hi', hi, 'open_left', i == 1, ...
                                     'open_right', j == n_g); %#ok<AGROW>
                i = j + 1;
            else
                i = i + 1;
            end
        end
        % Accepted intervals narrower than the grid step lie between two
        % rejected grid points: minimize the exact statistic around every
        % interior local minimum of the sampled one (lowest first, at most 50).
        if n_g >= 3
            ii = 2:n_g - 1;
            is_min = (stat_grid(ii) <= stat_grid(ii - 1)) & (stat_grid(ii) <= stat_grid(ii + 1)) ...
                     & ~accept(ii) & isfinite(stat_grid(ii));
            cand = ii(is_min);
            [~, ordr] = sort(stat_grid(cand));      % stable sort
            cand = cand(ordr);
            cand = cand(1:min(50, numel(cand)));
            for i_ = cand
                a_ = grid(i_ - 1);
                b_ = grid(i_ + 1);
                xatol = 1e-12 * max([abs(a_), abs(b_), sc]);
                [xm, fm] = fminbnd(f_, a_, b_, optimset('TolX', xatol, 'MaxIter', 500, ...
                                                        'MaxFunEvals', 1000, 'Display', 'off'));
                if isfinite(fm) && fm < 0
                    lo_ = ivm_ar_root(f_, a_, xm, sc);
                    hi_ = ivm_ar_root(f_, xm, b_, sc);
                    if ~isempty(lo_) && ~isempty(hi_)
                        segs(end+1) = struct('lo', lo_, 'hi', hi_, 'open_left', false, ...
                                             'open_right', false); %#ok<AGROW>
                    end
                end
            end
            if numel(segs) > 1
                [~, ordr] = sort([segs.lo]);
                segs = segs(ordr);
            end
        end
        out.segments = segs;
        if isempty(segs)
            out.shape = 'empty';
        else
            out.shape = 'grid';
        end
    end
end

% -- Does the set contain zero? --
switch out.shape
    case 'unavailable'
        out.contains_zero = [];
        return;
    case 'interval'
        out.contains_zero = (out.lower <= 0.0) && (0.0 <= out.upper);
    case 'rays'
        out.contains_zero = (0.0 <= out.left_upper) || (0.0 >= out.right_lower);
    case 'rays_left_only'
        out.contains_zero = (0.0 <= out.left_upper);
    case 'rays_right_only'
        out.contains_zero = (0.0 >= out.right_lower);
    case 'line'
        out.contains_zero = true;
    case 'grid'
        cz = false;
        for i = 1:numel(out.segments)
            if out.segments(i).lo <= 0.0 && 0.0 <= out.segments(i).hi
                cz = true;
                break;
            end
        end
        out.contains_zero = cz;
    otherwise   % 'empty'
        out.contains_zero = false;
end
% Consistency: the AR test at zero must agree with set membership (up to
% boundary ties); prefer the direct test when they disagree.
if isfinite(stat0)
    test_says_in = (stat0 <= cv);
    if test_says_in ~= out.contains_zero
        out.contains_zero = test_says_in;
    end
end
end

function [s, ev] = ivm_ar_closed(cf, b0)
% Closed-form AR statistic b(b0)' V(b0)^{-1} b(b0) (and the eigenvalues of
% the symmetrized V(b0)).
s_ = b0 - cf.tc;
bb = cf.b_c - s_ * cf.Bv;
VV = cf.V_c + s_ * cf.V1 + s_ * s_ * cf.V2;
if nargout > 1
    ev = eig((VV + VV') / 2);
end
if any(~isfinite(VV(:))) || any(~isfinite(bb))
    s = NaN;
    return;
end
[Lv, Uv, Pv] = lu(VV);
if any(diag(Uv) == 0)
    s = NaN;
    return;
end
s = bb' * (Uv \ (Lv \ (Pv * bb)));
end

function r = ivm_ar_root(f_, a_, b_, sc)
% Brent root of f_ on [a_, b_] (empty if no sign change).
fa = f_(a_);
fb = f_(b_);
if ~(isfinite(fa) && isfinite(fb)) || fa * fb > 0
    r = [];
    return;
end
if fa == 0
    r = a_;
    return;
end
if fb == 0
    r = b_;
    return;
end
xtol = 1e-14 * max([abs(a_), abs(b_), sc]);
r = fzero(f_, [a_, b_], optimset('TolX', xtol, 'Display', 'off'));
end

function [bz, Vz] = ivm_ar_eval(b0, y_dep, x_endog, W_incl, z_cols, sekw, df_adjust)
u = y_dep - b0 * x_endog;
res = ivm_twosls_with_vcov(u, W_incl, [], sekw, df_adjust);
bz = res.coef(z_cols);
Vz = res.vcov(z_cols, z_cols);
end

function s = ivm_ar_stat(b0, y_dep, x_endog, W_incl, z_cols, sekw, df_adjust)
[bz, Vz] = ivm_ar_eval(b0, y_dep, x_endog, W_incl, z_cols, sekw, df_adjust);
s = bz' * (Vz \ bz);
end

function tf = ivm_ar_unbounded_or_disconnected(ar)
% True when an AR set is unbounded (rays, line, or a grid segment reaching
% the grid boundary) or has more than one segment.
tf = false;
if isempty(ar)
    return;
end
switch ar.shape
    case {'rays', 'rays_left_only', 'rays_right_only', 'line'}
        tf = true;
    case 'grid'
        segs = ar.segments;
        tf = numel(segs) > 1 || any([segs.open_left]) || any([segs.open_right]);
end
end

function s = ivm_format_ar_set(ar, fmtstr)
%IVM_FORMAT_AR_SET Human-readable rendering of an AR confidence set.
if nargin < 2
    fmtstr = '%.4f';
end
f = @(v) sprintf(fmtstr, v);
switch ar.shape
    case 'interval'
        s = sprintf('[%s, %s]', f(ar.lower), f(ar.upper));
    case 'rays'
        s = sprintf(['(-inf, %s] U [%s, +inf)   <- unbounded: ' ...
                     'weak-IV signature'], f(ar.left_upper), f(ar.right_lower));
    case 'rays_left_only'
        s = sprintf('(-inf, %s]   <- unbounded: weak-IV signature', f(ar.left_upper));
    case 'rays_right_only'
        s = sprintf('[%s, +inf)   <- unbounded: weak-IV signature', f(ar.right_lower));
    case 'line'
        s = ['(-inf, +inf)   <- AR cannot reject any value at this level ' ...
             '(weak-IV signature)'];
    case 'grid'
        parts = cell(1, numel(ar.segments));
        any_open = false;
        for i = 1:numel(ar.segments)
            sg = ar.segments(i);
            if sg.open_left
                lo_s = ['< ', f(sg.lo)];
            else
                lo_s = f(sg.lo);
            end
            if sg.open_right
                hi_s = ['> ', f(sg.hi)];
            else
                hi_s = f(sg.hi);
            end
            any_open = any_open || sg.open_left || sg.open_right;
            parts{i} = sprintf('[%s, %s]', lo_s, hi_s);
        end
        if any_open
            note = ['   (grid search; ''<'' / ''>'' mark segments that reach the ' ...
                    'search grid''s boundary)'];
        else
            note = '';
        end
        s = [strjoin(parts, ' U '), note];
    case 'unavailable'
        s = 'unavailable: negative variance estimate';
    otherwise
        s = 'empty set (AR rejects every value: model misspecification signal)';
end
end

% ============================================================================
%  Variable lists (globs, Python fnmatch semantics, case-sensitive)
% ============================================================================

function rx = ivm_glob_to_regex(pat)
%IVM_GLOB_TO_REGEX Translate a shell glob to an anchored regex (fnmatch):
% '*' -> '.*', '?' -> '.', '[seq]' (with [!seq] negation, ranges, and a
% literal leading '^'); an unterminated '[' is a literal bracket.
chars = pat;
n = numel(chars);
i = 1;
out = '';
while i <= n
    ch = chars(i);
    if ch == '*'
        out = [out, '.*']; %#ok<AGROW>
    elseif ch == '?'
        out = [out, '.']; %#ok<AGROW>
    elseif ch == '['
        j = i + 1;
        if j <= n && chars(j) == '!'
            j = j + 1;
        end
        if j <= n && chars(j) == ']'
            j = j + 1;
        end
        while j <= n && chars(j) ~= ']'
            j = j + 1;
        end
        if j > n
            out = [out, '\[']; %#ok<AGROW>
        else
            inner = chars(i+1:j-1);
            if ~isempty(inner) && inner(1) == '!'
                inner = ['^', inner(2:end)];
            elseif ~isempty(inner) && inner(1) == '^'
                inner = ['\', inner];
            end
            out = [out, '[', inner, ']']; %#ok<AGROW>
            i = j;
        end
    else
        if any(ch == '][{}()*+?.\^$|')
            out = [out, '\', ch]; %#ok<AGROW>
        else
            out = [out, ch]; %#ok<AGROW>
        end
    end
    i = i + 1;
end
rx = ['^', out, '$'];
end

function tf = ivm_has_glob_chars(tok)
tf = any(tok == '*') || any(tok == '?') || any(tok == '[');
end

function hits = ivm_glob_filter(names, pat)
rx = ivm_glob_to_regex(pat);
keep = false(1, numel(names));
for i = 1:numel(names)
    keep(i) = ~isempty(regexp(names{i}, rx, 'once'));
end
hits = names(keep);
end

function hits = ivm_glob_filter_names(names, pat)
% Exported testing surface: glob-filter a cellstr of names.
hits = ivm_glob_filter(names, pat);
end

function tokens = ivm_spec_tokens(spec, name)
% Tokens of a variable spec: [] / '' -> {}; char or string scalar -> split
% on whitespace; cellstr / string array -> elements.
if isempty(spec) && ~ischar(spec)
    tokens = {};
elseif ischar(spec) || (isa(spec, 'string') && isscalar(spec))
    s = strtrim(ivm_norm_text_scalar(spec, name));
    if isempty(s)
        tokens = {};
    else
        tokens = regexp(s, '\s+', 'split');
    end
else
    tokens = ivm_norm_text_list(spec, name);
end
end

function [out, excluded] = ivm_resolve_var_list(spec, columns, data, name, exclude, forbid)
%IVM_RESOLVE_VAR_LIST Normalize a variable spec to a list of column names.
% Glob tokens (* ? [) expand against `columns`; literal tokens must exist.
% Columns in `exclude` (Y, T, M, Z and the fe/cluster/weight/coordinate
% columns) are never produced by a glob; a literal token naming a column in
% `forbid` (Y, T, M, Z) is an error. Duplicates are removed (first kept).
if nargin < 5, exclude = {}; end
if nargin < 6, forbid = {}; end
excluded = {};
tokens = ivm_spec_tokens(spec, name);
resolved = {};
for t = 1:numel(tokens)
    tok = tokens{t};
    if ivm_has_glob_chars(tok)
        hits = ivm_glob_filter(columns, tok);
        if isempty(hits)
            error(['%s: glob pattern ''%s'' matched 0 columns. ' ...
                   'Available columns include: %s%s'], name, tok, ...
                  strjoin(columns(1:min(8, numel(columns))), ', '), ...
                  ivm_ellipsis(numel(columns) > 8));
        end
        is_ex = ismember(hits, exclude);
        for m = find(is_ex)
            if ~any(strcmp(excluded, hits{m}))
                excluded{end+1} = hits{m}; %#ok<AGROW>
            end
        end
        if all(is_ex)
            error(['%s: glob pattern ''%s'' matched only role columns (%s), ' ...
                   'which cannot be controls.'], name, tok, strjoin(hits, ', '));
        end
        resolved = [resolved, hits(~is_ex)]; %#ok<AGROW>
    else
        if ~any(strcmp(columns, tok))
            error(['%s: variable ''%s'' not found in data. ' ...
                   'Available columns include: %s%s'], name, tok, ...
                  strjoin(columns(1:min(8, numel(columns))), ', '), ...
                  ivm_ellipsis(numel(columns) > 8));
        end
        resolved{end+1} = tok; %#ok<AGROW>
    end
end
bad = resolved(ismember(resolved, forbid));
if ~isempty(bad)
    error('%s must not contain the outcome, treatment, mediator or instruments: %s.', ...
          name, strjoin(ivm_unique_stable(bad), ', '));
end
out = ivm_unique_stable(resolved);
all_nan = {};
for i = 1:numel(out)
    col = ivm_df_get(ivm_df_from_user(data, out(i)), out{i});
    if all(ivm_col_isna(col))
        all_nan{end+1} = out{i}; %#ok<AGROW>
    end
end
if ~isempty(all_nan)
    error(['%s: column(s) %s contain only missing values. Either drop them ' ...
           'from %s or fill the missing values before calling ' ...
           'ivmediate.'], name, strjoin(all_nan, ', '), name);
end
end

function u = ivm_unique_stable(c)
u = {};
for i = 1:numel(c)
    if ~any(strcmp(u, c{i}))
        u{end+1} = c{i}; %#ok<AGROW>
    end
end
end

function s = ivm_ellipsis(cond)
if cond
    s = '...';
else
    s = '';
end
end

function line = ivm_summarize_resolved(name, resolved, max_show)
if nargin < 3
    max_show = 8;
end
if isempty(resolved)
    line = '';
    return;
end
n = numel(resolved);
if n <= max_show
    body = strjoin(resolved, ', ');
else
    body = sprintf('%s, ... (+%d more)', strjoin(resolved(1:max_show), ', '), n - max_show);
end
if n ~= 1
    plural = 's';
else
    plural = '';
end
line = sprintf('  Resolved %s (%d var%s): %s', name, n, plural, body);
end


% ============================================================================
%  Roles, missing values and the estimation sample (Python _prepare_sample)
% ============================================================================

function Zl = ivm_normalize_Z(Z)
% Z as a name, cellstr or string array of distinct instrument names.
if ischar(Z) || (isa(Z, 'string') && isscalar(Z))
    Zl = {ivm_norm_text_scalar(Z, 'Z')};
elseif iscell(Z) || isa(Z, 'string')
    Zl = ivm_norm_text_list(Z, 'instrument names');
else
    error('Z must be a column name or a list of names; got %s.', ivm_repr(Z));
end
if isempty(Zl)
    error('Z must name at least one instrument.');
end
dup = {};
for i = 2:numel(Zl)
    if any(strcmp(Zl(1:i-1), Zl{i})) && ~any(strcmp(dup, Zl{i}))
        dup{end+1} = Zl{i}; %#ok<AGROW>
    end
end
if ~isempty(dup)
    error('duplicate instrument names: %s.', strjoin(dup, ', '));
end
end

function fe_cols = ivm_normalize_fe(fe)
if isempty(fe) && ~ischar(fe)
    fe_cols = {};
    return;
end
fe_cols = ivm_norm_text_list(fe, 'fe');
fe_cols = fe_cols(~cellfun('isempty', fe_cols));
if numel(fe_cols) > 2
    error('At most two-way fixed effects are supported.');
end
end

function prep = ivm_prepare_sample(data, Y, T, M, Z, X, dummies, fe, se, clusterid, weights, lat, lon)
%IVM_PREPARE_SAMPLE Resolve and validate the variable roles, normalize
% missing values and return the estimation sample (before singleton removal
% and FE absorption). Role checks; X/dummies never take a role column by
% glob, and an explicit Y/T/M/Z among them is an error; text id columns (fe,
% cluster, dummies): '' and whitespace-only are missing; listwise deletion
% over every used column (cluster ids only with se='cluster', coordinates
% only with se='spatial'); +-Inf in a numeric role is an error; negative
% weights are an error; zero-weight rows are dropped.
notes = {};
Y = ivm_norm_text_scalar(Y, 'Y');
T = ivm_norm_text_scalar(T, 'T');
M = ivm_norm_text_scalar(M, 'M');
Zl = ivm_normalize_Z(Z);
if numel(unique({Y, T, M})) < 3
    error('Y, T and M must be three different columns.');
end
clash = Zl(ismember(Zl, {Y, T, M}));
if ~isempty(clash)
    error('instruments must differ from the outcome, treatment and mediator: %s.', ...
          strjoin(clash, ', '));
end
fe_cols = ivm_normalize_fe(fe);
if strcmp(se, 'cluster')
    cl_eff = clusterid;
else
    cl_eff = [];
end
if strcmp(se, 'spatial') && (isempty(lat) || isempty(lon))
    error('se=''spatial'' requires lat and lon column names.');
end
roles = [{Y, T, M}, Zl];
names = ivm_data_names(data);
need = roles;
need = [need, fe_cols];
if ~isempty(cl_eff), need{end+1} = cl_eff; end
if ~isempty(weights), need{end+1} = weights; end
if strcmp(se, 'spatial'), need = [need, {lat, lon}]; end
missing_cols = ivm_unique_stable(need(~ismember(need, names)));
if ~isempty(missing_cols)
    error('column(s) not found in data: %s.', strjoin(missing_cols, ', '));
end
aux = fe_cols;
for c = {clusterid, weights, lat, lon}
    if ~isempty(c{1}), aux{end+1} = c{1}; end %#ok<AGROW>
end
exclude = [roles, aux];
[X_res, ex_x] = ivm_resolve_var_list(X, names, data, 'X', exclude, roles);
[D_res, ex_d] = ivm_resolve_var_list(dummies, names, data, 'dummies', exclude, roles);
if ~isempty(ex_x)
    notes{end+1} = sprintf('  Note: X globs also matched role columns, excluded: %s.', strjoin(ex_x, ', '));
end
if ~isempty(ex_d)
    notes{end+1} = sprintf('  Note: dummies globs also matched role columns, excluded: %s.', strjoin(ex_d, ', '));
end
numeric_roles = [roles, X_res];
if ~isempty(weights), numeric_roles{end+1} = weights; end
if strcmp(se, 'spatial'), numeric_roles = [numeric_roles, {lat, lon}]; end
numeric_roles = ivm_unique_stable(numeric_roles);
id_cols = fe_cols;
if ~isempty(cl_eff), id_cols{end+1} = cl_eff; end
id_cols = ivm_unique_stable([id_cols, D_res]);
used = ivm_unique_stable([numeric_roles, id_cols]);
df = ivm_df_from_user(data, used);
for i = 1:numel(id_cols)
    c = id_cols{i};
    col = ivm_df_get(df, c);
    if iscell(col)
        blank = cellfun(@(s) isempty(strtrim(s)), col);
        if any(blank)
            col(blank) = {''};
            df = ivm_df_set(df, c, col);
        end
    end
end
for i = 1:numel(numeric_roles)
    c = numeric_roles{i};
    col = ivm_df_get(df, c);
    if ~isnumeric(col)
        % examples: the values that are not numbers (Python lists the first
        % three values, which may all be numbers)
        vals = col(~cellfun('isempty', col));
        bad = vals(isnan(str2double(vals)));
        if isempty(bad)
            bad = vals;
        end
        ex = ivm_unique_stable(bad);
        ex = ex(1:min(3, numel(ex)));
        ex = cellfun(@(s) ['''', s, ''''], ex, 'UniformOutput', false);
        error('column ''%s'' must be numeric; found values such as %s.', c, strjoin(ex, ', '));
    end
end
n0 = ivm_df_nrows(df);
mask = true(n0, 1);
for i = 1:numel(used)
    mask = mask & ~ivm_col_isna(ivm_df_get(df, used{i}));
end
inf_cols = {};
for i = 1:numel(numeric_roles)
    col = ivm_df_get(df, numeric_roles{i});
    if any(isinf(col(mask)))
        inf_cols{end+1} = numeric_roles{i}; %#ok<AGROW>
    end
end
if ~isempty(inf_cols)
    error('infinite values in column(s): %s.', strjoin(inf_cols, ', '));
end
n_zero = 0;
if ~isempty(weights)
    wv = ivm_df_get(df, weights);
    if any(wv(mask) < 0)
        error('Weights must be non-negative.');
    end
    zero = mask & (wv == 0);
    n_zero = sum(zero);
    if n_zero > 0
        if n_zero == sum(mask)
            error(['All weights are zero - no observations have positive ' ...
                   'weight to estimate from.']);
        end
        mask = mask & ~zero;
    end
end
df = ivm_df_rows(df, mask);
if ivm_df_nrows(df) == 0
    error('No observations remaining after dropping missing values.');
end
prep = struct();
prep.df = df;
prep.mask = mask;
prep.Y = Y; prep.T = T; prep.M = M;
prep.Z = Zl;
prep.X = X_res;
prep.dummies = D_res;
prep.fe_cols = fe_cols;
prep.clusterid = cl_eff;
prep.notes = notes;
prep.n_zero_weights = n_zero;
prep.used = used;
end


% ============================================================================
%  Fixed effects: singletons, weighted demeaning, degrees of freedom
% ============================================================================

function [df, n_groups, n_dropped] = ivm_drop_singletons(df, fe_col)
% Iteratively drop singleton FE groups until none remain (reghdfe convention).
n_before = ivm_df_nrows(df);
while true
    [gidx, G] = ivm_factorize(ivm_df_get(df, fe_col));
    counts = accumarray(gidx, 1, [G, 1]);
    is_singleton = counts(gidx) == 1;
    if ~any(is_singleton)
        break;
    end
    df = ivm_df_rows(df, ~is_singleton);
end
[~, n_groups] = ivm_factorize(ivm_df_get(df, fe_col));
n_dropped = n_before - ivm_df_nrows(df);
end

function M = ivm_group_means(S, sw, codes, V)
% Weighted group means of the columns of V, expanded back to rows
% (S = sparse aggregation matrix of the weights, sw = group weight sums).
gm = (S * V) ./ sw;
M = full(gm(codes, :));
end

function df = ivm_demean_within_groups(df, cols, fe_col, w)
% Demean `cols` within groups of fe_col: x - sum_g(w x) / sum_g(w)
% (ordinary means when w is all ones).
[codes, G] = ivm_factorize(ivm_df_get(df, fe_col));
n = numel(codes);
S = sparse(codes, 1:n, w, G, n);
sw = full(S * ones(n, 1));
V = ivm_extract_mat(df, cols);
V = V - ivm_group_means(S, sw, codes, V);
for c = 1:numel(cols)
    df = ivm_df_set(df, cols{c}, V(:, c));
end
end

function [df, converged] = ivm_demean_twoway(df, cols, fe_cols, w, quiet, tol, max_iter, max_direct)
% Two-way FE demeaning (weighted; Python _demean_twoway). When the FE
% dimension with fewer levels has at most max_direct (3,000) levels, the
% exact projection of ivm_demean_twoway_direct (Addendum A.2, O-2:
% alternating projections can stall on weakly connected designs such as a
% worker-firm chain). Otherwise alternating (weighted) projections with
% Irons-Tuck acceleration; convergence: the largest change in one iteration
% divided by the SD of the original column is below tol (1e-12) for every
% column that varies; at most max_iter (10,000) iterations.
if nargin < 6 || isempty(tol), tol = 1e-12; end
if nargin < 7 || isempty(max_iter), max_iter = 10000; end
if nargin < 8 || isempty(max_direct), max_direct = 3000; end
[c1, G1] = ivm_factorize(ivm_df_get(df, fe_cols{1}));
[c2, G2] = ivm_factorize(ivm_df_get(df, fe_cols{2}));
n = numel(c1);
V = ivm_extract_mat(df, cols);
converged = true;
if size(V, 2) == 0
    return;
end
R = ivm_demean_twoway_direct(V, c1, G1, c2, G2, w, max_direct);
if ~isempty(R)
    for c = 1:numel(cols)
        df = ivm_df_set(df, cols{c}, R(:, c));
    end
    return;
end
S1 = sparse(c1, 1:n, w, G1, n);
S2 = sparse(c2, 1:n, w, G2, n);
sw1 = full(S1 * ones(n, 1));
sw2 = full(S2 * ones(n, 1));
sd = std(V, 1, 1);
active = sd > 1e-14 * max(max(abs(V), [], 1), 1e-300);
scale = max(sd, 1e-300);
sweep = @(A) ivm_twoway_sweep(A, S1, sw1, c1, S2, sw2, c2);
converged = ~any(active);
change = 0;
it = 0;
while ~converged && it < max_iter
    F = sweep(V);
    it = it + 1;
    dF = abs(F(:, active) - V(:, active));
    change = max(max(dF, [], 1) ./ scale(active));
    if change < tol
        V = F;
        converged = true;
        break;
    end
    Gs = sweep(F);
    it = it + 1;
    D1 = F - V;
    D2 = Gs - F;
    DD = D2 - D1;
    den = sum(DD .* DD, 1);
    coef = zeros(1, size(V, 2));
    pos = den > 0;
    coef(pos) = sum(D2(:, pos) .* DD(:, pos), 1) ./ den(pos);
    V = Gs - coef .* D2;                    % Irons-Tuck extrapolation
end
if ~converged
    ivm_note(quiet, sprintf(['  Warning: two-way FE demeaning did not converge after ' ...
        '%d iterations (last relative change = %.2e > tol = %.0e). Results may be ' ...
        'imprecise; check for nearly-collinear FE dimensions.'], max_iter, change, tol));
end
for c = 1:numel(cols)
    df = ivm_df_set(df, cols{c}, V(:, c));
end
end

function A = ivm_twoway_sweep(A, S1, sw1, c1, S2, sw2, c2)
A = A - ivm_group_means(S1, sw1, c1, A);
A = A - ivm_group_means(S2, sw2, c2, A);
end

function R = ivm_demean_twoway_direct(V, c1, G1, c2, G2, w, max_small)
%IVM_DEMEAN_TWOWAY_DIRECT Exact weighted two-way FE residualization
% (Addendum A.2, O-2; the algorithm of Stata's _ivm_demean2_direct and
% Python's _demean_twoway_direct). Model V = alpha_A + beta_B + R by weighted
% least squares, with A the FE dimension with more levels and B the other.
% Eliminating alpha gives L beta = r for the B effects, with the Schur
% complement L = diag(s_B) - N' diag(1/s_A) N (s_A, s_B: weight sums per
% level; N(a, b): weight of the rows in level a of A and level b of B;
% r: B-level sums of w (V - within-A means)). L has one null direction per
% connected component of the FE graph, so one B level per component (the
% first in code order) is fixed at 0; one Cholesky factorization serves all
% columns (min-norm least squares if it fails). The residual is
% R = (V - mean_A V) - (beta_B - mean_A beta_B). Returns [] when the smaller
% dimension has more than max_small levels.
if G2 > G1
    cA = c2; cB = c1; Ga = G2; Gb = G1;
else
    cA = c1; cB = c2; Ga = G1; Gb = G2;
end
R = [];
if Gb > max_small
    return;
end
n = numel(cA);
w = w(:);
SA = sparse(cA, 1:n, w, Ga, n);
sA = full(SA * ones(n, 1));
sB = accumarray(cB, w, [Gb, 1]);
N = sparse(cA, cB, w, Ga, Gb);                          % duplicates are summed
Ns = spdiags(1 ./ sqrt(sA), 0, Ga, Ga) * N;
L = diag(sB) - full(Ns' * Ns);
L = 0.5 * (L + L');
lab = ivm_component_labels(Ga + Gb, cA, Ga + cB);
[~, first] = unique(lab(Ga + 1:end), 'first');          % first B level per component
keep = true(Gb, 1);
keep(first) = false;
V1 = V - ivm_group_means(SA, sA, cA, V);
rhs = full(sparse(cB, 1:n, w, Gb, n) * V1);
beta = zeros(Gb, size(V, 2));
kb = find(keep);
if ~isempty(kb)
    Lk = L(kb, kb);
    [U, flag] = chol(Lk);
    if flag == 0
        beta(kb, :) = U \ (U' \ rhs(kb, :));
    else
        beta(kb, :) = ivm_lstsq(Lk, rhs(kb, :));
    end
end
Bv = beta(cB, :);
R = V1 - (Bv - ivm_group_means(SA, sA, cA, Bv));
end

function [d, d_cr1, n_groups, nested] = ivm_fe_df(df, fe_cols, cluster_col)
%IVM_FE_DF Absorbed FE degrees of freedom beyond the intercept (rank of the
% FE dummies): one-way G - 1; two-way G1 + G2 - c - 1 with c the connected
% components of the bipartite FE graph. With cluster SEs, an FE dimension
% nested in the clusters (every group within one cluster) does not count
% toward CR1 (reghdfe/fixest); with one of two dimensions nested, the other
% counts G_other - 1 (reghdfe, pyfixest default).
nd = numel(fe_cols);
codes = cell(1, nd);
G = zeros(1, nd);
for i = 1:nd
    [codes{i}, G(i)] = ivm_factorize(ivm_df_get(df, fe_cols{i}));
end
nested = false(1, nd);
if ~isempty(cluster_col)
    clc = ivm_factorize(ivm_df_get(df, cluster_col));
    for i = 1:nd
        pairs = unique([codes{i}, clc], 'rows');
        nested(i) = size(pairs, 1) == numel(unique(pairs(:, 1)));
    end
end
if nd == 1
    d = G(1) - 1;
    if nested(1)
        d_cr1 = 0;
    else
        d_cr1 = d;
    end
else
    edges = unique([codes{1}, codes{2}], 'rows');
    c = ivm_n_components(G(1) + G(2), edges(:, 1), G(1) + edges(:, 2));
    d = G(1) + G(2) - c - 1;
    if nested(1) && nested(2)
        d_cr1 = 0;
    elseif nested(1)
        d_cr1 = G(2) - 1;
    elseif nested(2)
        d_cr1 = G(1) - 1;
    else
        d_cr1 = d;
    end
end
n_groups = sum(G);
end

function c = ivm_n_components(N, a, b)
% Number of connected components of an undirected graph with N nodes and
% edges (a(k), b(k)).
c = numel(unique(ivm_component_labels(N, a, b)));
end

function parent = ivm_component_labels(N, a, b)
% Component label (root node) of each of N nodes of an undirected graph with
% edges (a(k), b(k)), by union-find with path halving; duplicate edges are
% removed first.
e = unique([a(:), b(:)], 'rows');
a = e(:, 1);
b = e(:, 2);
parent = 1:N;
for k = 1:numel(a)
    x = a(k);
    while parent(x) ~= x
        parent(x) = parent(parent(x));
        x = parent(x);
    end
    y = b(k);
    while parent(y) ~= y
        parent(y) = parent(parent(y));
        y = parent(y);
    end
    if x ~= y
        parent(max(x, y)) = min(x, y);
    end
end
for k = 1:N
    x = k;
    while parent(x) ~= x
        x = parent(x);
    end
    parent(k) = x;
end
end


% ============================================================================
%  Dummy expansion (drop-first encoding)
% ============================================================================

function lab = ivm_cat_label(v)
% Label of a numeric category: the integer form for integral values, else the
% shortest '%.15g'/'%.16g'/'%.17g' form that reads back exactly (so distinct
% values get distinct labels; pandas prints float categories as '2.0').
if v == fix(v) && abs(v) < 1e17
    lab = sprintf('%.0f', v + 0);
    return;
end
for p = 15:17
    lab = sprintf(sprintf('%%.%dg', p), v);
    if str2double(lab) == v
        return;
    end
end
end

function [df, new_cols] = ivm_expand_dummies(df, dummy_vars)
%IVM_EXPAND_DUMMIES Indicators for all but the first (sorted) category;
% generated names must be unique and must not clash with existing columns.
new_cols = {};
for k = 1:numel(dummy_vars)
    var = dummy_vars{k};
    col = ivm_df_get(df, var);
    miss = ivm_col_isna(col);
    cats = sort(unique(col(~miss)));     % numeric sort, or code-point sort
    n_cats = numel(cats);
    if n_cats <= 1
        continue;
    end
    for c = 2:n_cats
        if isnumeric(col)
            cat_v = cats(c);
            ind = double(col == cat_v);
            lab = ivm_cat_label(cat_v);
        else
            cat_v = cats{c};
            ind = double(strcmp(col, cat_v));
            lab = cat_v;
        end
        col_name = ['_D_', var, '_', lab];
        if any(strcmp(new_cols, col_name)) || ivm_df_has(df, col_name)
            error(['dummies: the indicator name ''%s'' is not unique (another level ' ...
                   'or column has the same name); rename the levels or columns.'], col_name);
        end
        df = ivm_df_set(df, col_name, ind);
        new_cols{end+1} = col_name; %#ok<AGROW>
    end
end
end


% ============================================================================
%  Collinear columns (order-respecting, Python _greedy_independent)
% ============================================================================

function [kept, dropped] = ivm_greedy_independent(A, names, tol, keep_first, ref_norms)
% Columns of A that are linearly independent, tested in order: each column is
% divided by its reference scale and kept when the norm of its residual on
% the previously kept columns (two Gram-Schmidt passes) exceeds tol (1e-9).
% keep_first columns are orthogonalized first and always kept.
if nargin < 3 || isempty(tol), tol = 1e-9; end
[n, p] = size(A);
basis = zeros(n, 0);
if nargin >= 4 && ~isempty(keep_first)
    K0 = reshape(keep_first, n, []);
    for j = 1:size(K0, 2)
        v = K0(:, j) / max(norm(K0(:, j)), 1e-300);
        for b = 1:size(basis, 2)
            v = v - (basis(:, b)' * v) * basis(:, b);
        end
        nv = norm(v);
        if nv > tol
            basis(:, end+1) = v / nv; %#ok<AGROW>
        end
    end
end
kept = [];
dropped = {};
for j = 1:p
    if nargin >= 5 && ~isempty(ref_norms)
        ref = ref_norms(j);
    else
        ref = norm(A(:, j));
    end
    if ~isfinite(ref) || ref <= 0
        dropped{end+1} = names{j}; %#ok<AGROW>
        continue;
    end
    v = A(:, j) / ref;
    for pass = 1:2
        for b = 1:size(basis, 2)
            v = v - (basis(:, b)' * v) * basis(:, b);
        end
    end
    nv = norm(v);
    if isfinite(nv) && nv > tol
        basis(:, end+1) = v / nv; %#ok<AGROW>
        kept(end+1) = j; %#ok<AGROW>
    else
        dropped{end+1} = names{j}; %#ok<AGROW>
    end
end
end

function Zm = ivm_extract_mat(df, names)
n = ivm_df_nrows(df);
Zm = zeros(n, numel(names));
for i = 1:numel(names)
    Zm(:, i) = ivm_col_num(ivm_df_get(df, names{i}), names{i});
end
end

% ============================================================================
%  Rich expansion for post-double-selection LASSO
% ============================================================================

function tf = ivm_is_continuous(values)
% A variable is continuous if it has more than 2 unique values (indicators
% are not).
tf = numel(unique(values)) > 2;
end

function ops = ivm_ops_for(var, default_ops, override)
% Effective expansion ops for one variable, honoring expand_override.
if ~isempty(override) && override.isKey(var)
    ops = ivm_norm_text_list(override(var), 'expand_override');
else
    ops = default_ops;
end
end

function override = ivm_norm_override(expand_override)
% expand_override as a containers.Map (or empty): a containers.Map, an Nx2
% cell {name, {ops...}} or a struct (field names must be valid identifiers;
% use the Map or cell form for names such as '_D_dnum_2').
if isempty(expand_override)
    override = [];
    return;
end
if isa(expand_override, 'containers.Map')
    override = expand_override;
    return;
end
if iscell(expand_override) && size(expand_override, 2) == 2
    override = containers.Map('KeyType', 'char', 'ValueType', 'any');
    for i = 1:size(expand_override, 1)
        override(ivm_norm_text_scalar(expand_override{i, 1}, 'expand_override')) = ...
            expand_override{i, 2};
    end
    return;
end
if isstruct(expand_override) && isscalar(expand_override)
    override = containers.Map('KeyType', 'char', 'ValueType', 'any');
    fn = fieldnames(expand_override);
    for i = 1:numel(fn)
        override(fn{i}) = expand_override.(fn{i});
    end
    return;
end
error(['expand_override must be a containers.Map, an Nx2 cell ' ...
       '{name, {ops...}}, or a struct.']);
end

function [W, names, meta] = ivm_build_rich_expansion(df, x_names, ...
    dummy_col_names, expand_x, expand_dummies, expand_override, quiet)
%IVM_BUILD_RICH_EXPANSION Candidate set: originals + squares/cubes +
% interactions, per the user-controlled expansion operations (built on the
% level-form data, before FE demeaning and weighting).
VALID_OPS = {'square', 'cube', 'interact_with_X', 'interact_with_dummies'};
override = ivm_norm_override(expand_override);

all_ops = {expand_x, expand_dummies};
if ~isempty(override)
    ov_vals = override.values;
    for i = 1:numel(ov_vals)
        all_ops{end+1} = ivm_norm_text_list(ov_vals{i}, 'expand_override'); %#ok<AGROW>
    end
end
for i = 1:numel(all_ops)
    ops_i = all_ops{i};
    for j = 1:numel(ops_i)
        if ~any(strcmp(VALID_OPS, ops_i{j}))
            error('Unknown expansion op ''%s''. Valid: %s', ops_i{j}, ...
                  strjoin(sort(VALID_OPS), ', '));
        end
    end
end

n = ivm_df_nrows(df);
cont_x = {};
for i = 1:numel(x_names)
    if ivm_is_continuous(ivm_df_get(df, x_names{i}))
        cont_x{end+1} = x_names{i}; %#ok<AGROW>
    end
end

X_TIMES = [' ', ivm_uch_raw(215), ' '];   % ' x ' with U+00D7 (names are data)

cols = {};
names = {};
meta = struct();
meta.n_originals_x = numel(x_names);
meta.n_originals_dummies = numel(dummy_col_names);
meta.n_continuous_x = numel(cont_x);
meta.n_dummies = numel(dummy_col_names);
meta.n_squares = 0;
meta.n_cubes = 0;
meta.n_interactions_xx = 0;
meta.n_interactions_dx = 0;
meta.n_interactions_dd = 0;
meta.warnings = {};

% -- Originals (X variables + dummy indicators) come first --
for i = 1:numel(x_names)
    cols{end+1} = ivm_col_num(ivm_df_get(df, x_names{i}), x_names{i}); %#ok<AGROW>
    names{end+1} = x_names{i}; %#ok<AGROW>
end
for i = 1:numel(dummy_col_names)
    cols{end+1} = ivm_col_num(ivm_df_get(df, dummy_col_names{i}), dummy_col_names{i}); %#ok<AGROW>
    names{end+1} = dummy_col_names{i}; %#ok<AGROW>
end

% -- Squares and cubes (continuous X only) --
for i = 1:numel(x_names)
    v = x_names{i};
    ops = ivm_ops_for(v, expand_x, override);
    is_cont = any(strcmp(cont_x, v));
    if any(strcmp(ops, 'square'))
        if ~is_cont
            meta.warnings{end+1} = sprintf(['''square'' has no effect on ' ...
                'indicator-like variable ''%s''; skipped.'], v);
        else
            cols{end+1} = ivm_col_num(ivm_df_get(df, v), v).^2; %#ok<AGROW>
            names{end+1} = [v, '^2']; %#ok<AGROW>
            meta.n_squares = meta.n_squares + 1;
        end
    end
    if any(strcmp(ops, 'cube'))
        if ~is_cont
            meta.warnings{end+1} = sprintf(['''cube'' has no effect on ' ...
                'indicator-like variable ''%s''; skipped.'], v);
        else
            cols{end+1} = ivm_col_num(ivm_df_get(df, v), v).^3; %#ok<AGROW>
            names{end+1} = [v, '^3']; %#ok<AGROW>
            meta.n_cubes = meta.n_cubes + 1;
        end
    end
end

% -- Continuous-continuous interactions (pairs from x_names) --
for i = 1:numel(x_names)
    v = x_names{i};
    ops = ivm_ops_for(v, expand_x, override);
    if ~any(strcmp(ops, 'interact_with_X'))
        continue;
    end
    for j = i+1:numel(x_names)    % only (v, u) with v before u in x_names
        u = x_names{j};
        cols{end+1} = ivm_col_num(ivm_df_get(df, v), v) .* ivm_col_num(ivm_df_get(df, u), u); %#ok<AGROW>
        names{end+1} = [v, X_TIMES, u]; %#ok<AGROW>
        meta.n_interactions_xx = meta.n_interactions_xx + 1;
    end
end

% -- Dummy x X interactions --
for i = 1:numel(dummy_col_names)
    d = dummy_col_names{i};
    ops = ivm_ops_for(d, expand_dummies, override);
    if any(strcmp(ops, 'square')) || any(strcmp(ops, 'cube'))
        meta.warnings{end+1} = sprintf(['''square''/''cube'' on dummy ' ...
            'variable ''%s'' has no effect; skipped.'], d);
    end
    if any(strcmp(ops, 'interact_with_X'))
        for j = 1:numel(x_names)
            u = x_names{j};
            cols{end+1} = ivm_col_num(ivm_df_get(df, d), d) .* ivm_col_num(ivm_df_get(df, u), u); %#ok<AGROW>
            names{end+1} = [d, X_TIMES, u]; %#ok<AGROW>
            meta.n_interactions_dx = meta.n_interactions_dx + 1;
        end
    end
end

% -- Dummy x dummy interactions (skip same-source indicator pairs) --
for i = 1:numel(dummy_col_names)
    d = dummy_col_names{i};
    ops = ivm_ops_for(d, expand_dummies, override);
    if ~any(strcmp(ops, 'interact_with_dummies'))
        continue;
    end
    for j = i+1:numel(dummy_col_names)
        e = dummy_col_names{j};
        if strncmp(d, '_D_', 3) && strncmp(e, '_D_', 3)
            if strcmp(ivm_rsplit_head(d), ivm_rsplit_head(e))
                continue;   % same source categorical -> product is 0
            end
        end
        cols{end+1} = ivm_col_num(ivm_df_get(df, d), d) .* ivm_col_num(ivm_df_get(df, e), e); %#ok<AGROW>
        names{end+1} = [d, X_TIMES, e]; %#ok<AGROW>
        meta.n_interactions_dd = meta.n_interactions_dd + 1;
    end
end

if isempty(cols)
    W = zeros(n, 0);
else
    W = zeros(n, numel(cols));
    for i = 1:numel(cols)
        W(:, i) = cols{i};
    end
end

% Duplicate names are ambiguous for selection-by-name: raise.
sorted_names = sort(names);
dups = {};
for i = 2:numel(sorted_names)
    if strcmp(sorted_names{i}, sorted_names{i-1})
        if isempty(dups) || ~strcmp(dups{end}, sorted_names{i})
            dups{end+1} = sorted_names{i}; %#ok<AGROW>
        end
    end
end
if ~isempty(dups)
    error(['Rich expansion produced duplicate column names: %s. This ' ...
           'usually means a user-specified variable was already an ' ...
           'expansion-style term (e.g., ''age^2'') and the same expansion ' ...
           'would generate it again. Either remove the duplicate from X, ' ...
           'exclude the operation via expand_override, or rename the ' ...
           'user-supplied column.'], strjoin(dups, ', '));
end

% Drop near-constant columns (variance below tolerance).
if size(W, 2) > 0
    col_var = var(W, 1, 1);                          % ddof = 0, like np.var
    col_mean_abs = max(abs(mean(W, 1)), 1.0);
    keep_mask = col_var > (1e-12 * col_mean_abs.^2);
    n_dropped = sum(~keep_mask);
    if n_dropped > 0
        meta.warnings{end+1} = sprintf(['Dropped %d near-constant columns ' ...
            'from rich expansion.'], n_dropped);
    end
    W = W(:, keep_mask);
    names = names(keep_mask);
end

for i = 1:min(3, numel(meta.warnings))
    ivm_note(quiet, sprintf('  Note: %s', meta.warnings{i}));
end
if numel(meta.warnings) > 3
    ivm_note(quiet, sprintf('  Note: ... (+%d more expansion warnings)', ...
                            numel(meta.warnings) - 3));
end
meta.p = size(W, 2);
end

function head = ivm_rsplit_head(s)
% Everything before the LAST underscore (Python s.rsplit('_', 1)[0]).
k = find(s == '_', 1, 'last');
if isempty(k)
    head = s;
else
    head = s(1:k-1);
end
end


% ============================================================================
%  BCH plug-in LASSO (penalty loadings, effective n, coordinate descent)
%  Since 2026-10: weights normalized to mean 1, y and W centered at their
%  weighted means (unpenalized intercept); heteroskedastic loadings
%  psi_j = sqrt(sum (w W_j u)^2 / n_w); n_eff = n or Kish n in every regime;
%  post-LASSO residuals for the loadings (hdm/lassopack default).
% ============================================================================

function psi = ivm_compute_psi(W, u, weights, cl_gidx, cl_G, spatial_A)
%IVM_COMPUTE_PSI Per-feature penalty loadings from the current residuals;
% the score of observation i is w_i W_ij u_i.
n_w = sum(weights);
if ~isempty(spatial_A)
    wuW = W .* (weights .* u);
    psi2 = sum(wuW .* (spatial_A * wuW), 1).';
    psi = sqrt(max(psi2, 0.0) / n_w);
elseif isempty(cl_gidx)
    wu2 = (weights .* u) .^ 2;
    psi = sqrt(sum(W .* W .* wu2, 1).' / n_w);
else
    wuW = W .* (weights .* u);
    S = sparse(cl_gidx, 1:numel(u), 1, cl_G, numel(u));
    cluster_sums = full(S * wuW);
    psi = sqrt(sum(cluster_sums .^ 2, 1).' / n_w);
end
psi = max(psi, 1e-10 * max([1.0; abs(W(:))]));
end

function n_eff = ivm_compute_n_eff(weights)
%IVM_COMPUTE_N_EFF n, or the Kish size (sum w)^2 / sum w^2 with unequal
% weights, in every regime (dependence enters through the loadings only).
if ivm_all_same_weight(weights)
    n_eff = numel(weights);
else
    n_w = sum(weights);
    n_eff = n_w * n_w / sum(weights .^ 2);
end
end

function lam = ivm_bch_lambda(p, n_eff, c_slack, alpha_p)
% BCH plug-in penalty level: lambda = 2c sqrt(n_eff) Phi^-1(1 - alpha_p/(2p)).
if p == 0
    lam = 0.0;
    return;
end
lam = 2.0 * c_slack * sqrt(n_eff) * ivm_norminv(1.0 - alpha_p / (2.0 * p));
end

function beta = ivm_lasso_cd(X, y, alpha, max_iter, tol)
%IVM_LASSO_CD Cyclic coordinate descent replicating sklearn's
% enet_coordinate_descent for Lasso (l1_ratio=1, fit_intercept=False):
% minimize (1/(2n))||y - Xb||^2 + alpha*||b||_1. Includes sklearn's
% duality-gap stopping rule (gap tolerance scaled by ||y||^2).
if nargin < 4, max_iter = 20000; end
if nargin < 5, tol = 1e-6; end
[n, p] = size(X);
beta = zeros(p, 1);
norm_cols = sum(X .* X, 1);
R = y;                              % residual at beta = 0
alpha_n = alpha * n;
tol_gap = tol * (y' * y);
for it = 1:max_iter
    w_max = 0.0;
    d_w_max = 0.0;
    for j = 1:p
        nj = norm_cols(j);
        if nj == 0
            continue;
        end
        w_old = beta(j);
        tmp = X(:, j)' * R + nj * w_old;
        bj = sign(tmp) * max(abs(tmp) - alpha_n, 0) / nj;
        if bj ~= w_old
            R = R - X(:, j) * (bj - w_old);
        end
        beta(j) = bj;
        dw = abs(bj - w_old);
        if dw > d_w_max
            d_w_max = dw;
        end
        if abs(bj) > w_max
            w_max = abs(bj);
        end
    end
    if w_max == 0.0 || d_w_max / w_max < tol || it == max_iter
        XtR = X' * R;
        dual_norm_XtA = max(abs(XtR));
        R_norm2 = R' * R;
        if dual_norm_XtA > alpha_n
            const = alpha_n / dual_norm_XtA;
            gap = 0.5 * (R_norm2 + R_norm2 * const^2);
        else
            const = 1.0;
            gap = R_norm2;
        end
        gap = gap + alpha_n * sum(abs(beta)) - const * (R' * y);
        if gap < tol_gap
            break;
        end
    end
end
end

function coef = ivm_lstsq(A, B)
% Minimum-norm least squares (numpy.linalg.lstsq with rcond=None: singular
% values below eps*max(m,n)*s_max are treated as zero).
if isempty(A)
    coef = zeros(size(A, 2), size(B, 2));
    return;
end
[U, S, V] = svd(A, 'econ');
s = diag(S);
tol = eps * max(size(A)) * max([s; 0]);
r = sum(s > tol);
coef = V(:, 1:r) * ((U(:, 1:r)' * B) ./ s(1:r));
end

function r = ivm_post_lasso_resid(y_c, W_c, sel, sqrt_w)
% Residuals of a weighted OLS of y on the selected columns of W (y and W are
% centered with the keep columns partialed out, so this is the post-LASSO
% regression on [selected, keep, intercept] by FWL).
if ~any(sel)
    r = y_c;
    return;
end
coef = ivm_lstsq(W_c(:, sel) .* sqrt_w, y_c .* sqrt_w);
r = y_c - W_c(:, sel) * coef;
end

function out = ivm_bch_lasso(y_raw, W_raw, weights, cl_gidx, cl_G, opts)
%IVM_BCH_LASSO Iterative plug-in LASSO with weight/cluster/spatial support
% (Python _bch_lasso). opts fields: n_iter_max, tol, c_slack, alpha_p,
% lambda_method, keep_idx, lambda_override, spatial_A, verbose.
[n, p] = size(W_raw);
if isempty(opts.alpha_p)
    alpha_p = 0.1 / max(log(max(n, 2)), 1.0);
else
    alpha_p = opts.alpha_p;
end
if p == 0
    out = ivm_lasso_out(zeros(0, 1), false(0, 1), 0.0, ones(0, 1), 0, true, 'trivial');
    return;
end
if strcmp(opts.lambda_method, 'cv') && isempty(opts.lambda_override)
    error(['lasso_penalty=''cv'' is not available in the MATLAB port: ' ...
           'sklearn''s KFold/GroupKFold fold assignment cannot be ' ...
           'reproduced here, so the CV path cannot be validated against ' ...
           'the Python output. Use lasso_penalty=''plugin'' (the default, ' ...
           'and the paper''s choice), or pass lasso_lambda_override to fix ' ...
           'lambda directly.']);
end

% Weight scale and location do not matter (2026-10).
weights = weights(:) / mean(weights);
y_raw = y_raw(:) - sum(weights .* y_raw(:)) / sum(weights);
W_raw = W_raw - sum(W_raw .* weights, 1) / sum(weights);
sqrt_w = sqrt(weights);

keep_idx = opts.keep_idx(:).';
if ~isempty(keep_idx)
    nonkeep_idx = setdiff(1:p, keep_idx);
    K_pm = W_raw(:, keep_idx) .* sqrt_w;
    coef_y = ivm_lstsq(K_pm, y_raw .* sqrt_w);
    y_lasso_raw = y_raw - W_raw(:, keep_idx) * coef_y;
    W_nonkeep = W_raw(:, nonkeep_idx);
    coef_W = ivm_lstsq(K_pm, W_nonkeep .* sqrt_w);
    W_lasso_raw = W_nonkeep - W_raw(:, keep_idx) * coef_W;
else
    nonkeep_idx = 1:p;
    y_lasso_raw = y_raw;
    W_lasso_raw = W_raw;
end
p_lasso = size(W_lasso_raw, 2);
y_lasso_pm = y_lasso_raw .* sqrt_w;
W_lasso_pm = W_lasso_raw .* sqrt_w;

n_eff = ivm_compute_n_eff(weights);

if p_lasso == 0
    selected = false(p, 1);
    if ~isempty(keep_idx)
        selected(keep_idx) = true;
    end
    out = ivm_lasso_out(zeros(p, 1), selected, 0.0, ones(p, 1), 0, true, 'keep-only');
    return;
end

if ~isempty(opts.lambda_override)
    % Loadings once from the initial residuals (the centered y), then a single
    % LASSO at the given lambda.
    psi_lasso = ivm_compute_psi(W_lasso_raw, y_lasso_raw, weights, cl_gidx, cl_G, opts.spatial_A);
    col_scale = 1.0 ./ psi_lasso;
    lam = opts.lambda_override;
    beta_cd = ivm_lasso_cd(W_lasso_pm .* col_scale.', y_lasso_pm, lam / (2.0 * n));
    beta_unscaled = beta_cd .* col_scale;
    beta = zeros(p, 1);
    beta(nonkeep_idx) = beta_unscaled;
    psi_full = NaN(p, 1);
    psi_full(nonkeep_idx) = psi_lasso;
    selected = false(p, 1);
    selected(nonkeep_idx) = beta_unscaled ~= 0;
    if ~isempty(keep_idx)
        selected(keep_idx) = true;
    end
    out = ivm_lasso_out(beta, selected, lam, psi_full, 1, true, 'override');
    return;
end

% -- Iterative plug-in path --
u = y_lasso_raw;                       % initial residuals: centered y
beta_lasso = zeros(p_lasso, 1);
n_iter = 0;
converged = false;
lam = 0.0;
psi_lasso = ones(p_lasso, 1);
for it = 1:opts.n_iter_max
    n_iter = it;
    psi_lasso = ivm_compute_psi(W_lasso_raw, u, weights, cl_gidx, cl_G, opts.spatial_A);
    lam = ivm_bch_lambda(p_lasso, n_eff, opts.c_slack, alpha_p);
    col_scale = 1.0 ./ psi_lasso;
    beta_cd = ivm_lasso_cd(W_lasso_pm .* col_scale.', y_lasso_pm, lam / (2.0 * n));
    beta_new = beta_cd .* col_scale;
    u_new = ivm_post_lasso_resid(y_lasso_raw, W_lasso_raw, beta_new ~= 0, sqrt_w);
    max_change = max(abs(beta_new - beta_lasso));
    if opts.verbose
        ivm_out('    iter %d: %s=%.3f, |%sb|_inf=%.2e, active=%d/%d\n', it, ...
                ivm_uch(955), lam, ivm_uch(916), max_change, sum(beta_new ~= 0), p_lasso);
    end
    beta_lasso = beta_new;
    u = u_new;
    if max_change < opts.tol
        converged = true;
        break;
    end
end

beta = zeros(p, 1);
beta(nonkeep_idx) = beta_lasso;
psi_full = NaN(p, 1);
psi_full(nonkeep_idx) = psi_lasso;
selected = false(p, 1);
selected(nonkeep_idx) = beta_lasso ~= 0;
if ~isempty(keep_idx)
    selected(keep_idx) = true;
end
out = ivm_lasso_out(beta, selected, lam, psi_full, n_iter, converged, 'plugin');
end

function out = ivm_lasso_out(beta, selected, lambda_value, psi, n_iter, converged, method)
% Scalar-struct builder (struct() with vector args would create struct arrays).
out = struct();
out.beta = beta;
out.selected = selected;
out.lambda_value = lambda_value;
out.psi = psi;
out.n_iter = n_iter;
out.converged = converged;
out.method = method;
end

function out = ivm_quadruple_lasso(y_raw, t_raw, m_raw, z_mat_raw, W_raw, ...
    col_names, weights, cl_gidx, cl_G, keep_idx, lopts)
%IVM_QUADRUPLE_LASSO Four LASSOs (Y, T, M, Z[k]) on the shared rich
% expansion; union of the selections plus per-LASSO diagnostics. M is
% regressed on W alone (paper App. H, Algorithm 1).
K_inst = size(z_mat_raw, 2);
p = size(W_raw, 2);
bopts = struct();
bopts.n_iter_max = lopts.n_iter_max;
bopts.tol = lopts.tol;
bopts.c_slack = lopts.c_slack;
bopts.alpha_p = lopts.alpha_p;
bopts.lambda_method = lopts.lambda_method;
bopts.keep_idx = keep_idx;
bopts.lambda_override = lopts.lambda_override;
bopts.spatial_A = lopts.spatial_A;
bopts.verbose = lopts.verbose;

if lopts.verbose
    ivm_out('  Running quadruple-LASSO:\n');
    ivm_out('    p (rich expansion) = %d\n', p);
    ivm_out('    K (instruments)    = %d\n', K_inst);
    ivm_out('    keep variables     = %d\n', numel(keep_idx));
end

if lopts.verbose, ivm_out('  %s LASSO Y on W\n', ivm_uch(8594)); end
res_Y = ivm_bch_lasso(y_raw, W_raw, weights, cl_gidx, cl_G, bopts);
if lopts.verbose, ivm_out('  %s LASSO T on W\n', ivm_uch(8594)); end
res_T = ivm_bch_lasso(t_raw, W_raw, weights, cl_gidx, cl_G, bopts);
if lopts.verbose
    ivm_out('  %s LASSO M on W (M ~ W only %s see Remark on omitting T)\n', ...
            ivm_uch(8594), ivm_uch(8212));
end
res_M = ivm_bch_lasso(m_raw, W_raw, weights, cl_gidx, cl_G, bopts);
res_Z_per_inst = cell(1, K_inst);
for k = 1:K_inst
    if lopts.verbose, ivm_out('  %s LASSO Z[%d] on W\n', ivm_uch(8594), k - 1); end
    res_Z_per_inst{k} = ivm_bch_lasso(z_mat_raw(:, k), W_raw, weights, ...
                                      cl_gidx, cl_G, bopts);
end

S_Y = res_Y.selected;
S_T = res_T.selected;
S_M = res_M.selected;
S_Z = false(p, 1);
for k = 1:K_inst
    S_Z = S_Z | res_Z_per_inst{k}.selected;
end
S = S_Y | S_T | S_M | S_Z;
if ~isempty(keep_idx)
    S(keep_idx) = true;
end

details = struct();
details.Y = ivm_lasso_summarize(res_Y, 'Y', col_names);
details.T = ivm_lasso_summarize(res_T, 'T', col_names);
details.M = ivm_lasso_summarize(res_M, 'M', col_names);
details.Z = cell(1, K_inst);
for k = 1:K_inst
    details.Z{k} = ivm_lasso_summarize(res_Z_per_inst{k}, ...
                                       sprintf('Z[%d]', k - 1), col_names);
end
S_Z_per_inst = cell(1, K_inst);
for k = 1:K_inst
    S_Z_per_inst{k} = res_Z_per_inst{k}.selected;
end

out = struct();
out.S = S;
out.S_Y = S_Y;
out.S_T = S_T;
out.S_M = S_M;
out.S_Z = S_Z;
out.S_Z_per_inst = S_Z_per_inst;
out.details = details;
out.union_names = col_names(S(:).');
out.n_selected_union = sum(S);
out.n_total = p;
end

function d = ivm_lasso_summarize(res, label, col_names)
sel = res.selected;
idx = find(sel(:).');
d = struct();
d.label = label;
d.names = col_names(idx);
d.coefs = res.beta(sel);
d.indices = idx;                        % 1-based positions in the expansion
d.lambda_value = res.lambda_value;
d.n_selected = sum(sel);
d.converged = res.converged;
d.n_iter = res.n_iter;
d.method = res.method;
end

function keep = ivm_resolve_lasso_keep(lasso_keep, rich_col_names)
%IVM_RESOLVE_LASSO_KEEP Resolve lasso_keep against the rich-expansion names
% (originals, powers 'x1^2' and interactions 'x1 <U+00D7> x2'; globs allowed).
tokens = ivm_spec_tokens(lasso_keep, 'lasso_keep');
keep = {};
for t = 1:numel(tokens)
    tok = tokens{t};
    if ivm_has_glob_chars(tok)
        hits = ivm_glob_filter(rich_col_names, tok);
        if isempty(hits)
            error(['lasso_keep: glob pattern ''%s'' matched 0 names in the ' ...
                   'rich expansion. Available names include: %s%s'], tok, ...
                  strjoin(rich_col_names(1:min(8, numel(rich_col_names))), ', '), ...
                  ivm_ellipsis(numel(rich_col_names) > 8));
        end
        for m = 1:numel(hits)
            if ~any(strcmp(keep, hits{m}))
                keep{end+1} = hits{m}; %#ok<AGROW>
            end
        end
    else
        if ~any(strcmp(rich_col_names, tok))
            error(['lasso_keep variable ''%s'' is not in the rich expansion. ' ...
                   'lasso_keep accepts originals (e.g., ''X1''), powers ' ...
                   '(e.g., ''X1^2''), and interactions (e.g., ''X1 %s X2''). ' ...
                   'Available names include: %s%s'], tok, ivm_uch(215), ...
                  strjoin(rich_col_names(1:min(8, numel(rich_col_names))), ', '), ...
                  ivm_ellipsis(numel(rich_col_names) > 8));
        end
        if ~any(strcmp(keep, tok))
            keep{end+1} = tok; %#ok<AGROW>
        end
    end
end
end

% ============================================================================
%  Main analysis function
% ============================================================================

function results = ivm_analysis(data, Y, T, M, Z, varargin)
%IVM_ANALYSIS Full IV mediation analysis with diagnostics (Python parity).
% Public entry point: output mode ('ascii') and MATLAB warnings are set for
% the duration of the call and restored afterward.
defaults = ivm_analysis_defaults();
opts = ivm_parse_opts(defaults, varargin, 'ivmediate');
[old_ascii, wstate] = ivm_enter(opts.ascii);
cleanup = onCleanup(@() ivm_leave(old_ascii, wstate)); %#ok<NASGU>
results = ivm_analysis_impl(data, Y, T, M, Z, opts, false);
end

function d = ivm_analysis_defaults()
d = struct();
d.X = [];
d.dummies = [];
d.fe = [];
d.se = 'robust';
d.clusterid = [];
d.alpha = 0.05;
d.weights = [];
d.wtype = 'aweight';
d.op_tau = 0.10;
d.op_alpha = 0.05;
d.ar_ci = true;
d.ar_grid_halfwidth = 200.0;
d.ar_grid_points = 4001;
d.lat = [];
d.lon = [];
d.spatial_cutoff = [];
d.spatial_kernel = [];        % 'uniform' with se='spatial'
d.spatial_metric = [];        % 'haversine' with se='spatial'
d.lasso = false;
d.lasso_penalty = 'plugin';
d.lasso_keep = [];
d.lasso_seed = 20260317;
d.lasso_cv_folds = 5;
d.expand_x = {'square', 'cube'};
d.expand_dummies = {'interact_with_X'};
d.expand_override = [];
d.show_selection = false;
d.lasso_c_slack = 1.1;
d.lasso_alpha_p = [];
d.lasso_n_iter_max = 15;
d.lasso_tol = 1e-6;
d.lasso_lambda_override = [];
d.lasso_verbose = false;
d.quiet = false;
d.ascii = [];                 % [] = automatic (on in matlab -batch on Windows)
end

function [old_ascii, wstate] = ivm_enter(ascii_opt)
% Set the output mode and silence MATLAB's singular-matrix warnings (the
% code checks conditioning itself; the warnings would only flood the log).
old_ascii = ivm_ascii_state();
new_ascii = ivm_resolve_ascii(ascii_opt);
wstate = warning();
warning('off', 'MATLAB:singularMatrix');
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:rankDeficientMatrix');
warning('off', 'MATLAB:illConditionedMatrix');
ivm_ascii_state(new_ascii);
end

function ivm_leave(old_ascii, wstate)
ivm_ascii_state(old_ascii);
warning(wstate);
end

function name = ivm_opt_name(v, what)
% Optional column-name option: [] (or '') stays [], text becomes char.
if nargin < 2
    what = 'option';
end
if (isempty(v) && ~isa(v, 'string')) || (isa(v, 'string') && isscalar(v) && strlength(v) == 0)
    name = [];
else
    name = ivm_norm_text_scalar(v, what);
    if isempty(name)
        name = [];
    end
end
end

function results = ivm_analysis_impl(data, Y, T, M, Z, opts, point_only)
%IVM_ANALYSIS_IMPL The analysis proper. point_only (internal, bootstrap
% replicates of basic/percentile intervals) skips every variance, the AR
% sets and the spatial neighbor lists (unless the LASSO needs them).
quiet = ivm_check_bool('quiet', opts.quiet);

% -- Validate options (2026-10: enumerations, ranges and booleans) --
se = ivm_check_choice('se', opts.se, {'classical', 'robust', 'cluster', 'spatial'});
wtype = ivm_check_choice('wtype', opts.wtype, {'aweight', 'pweight'});
ivm_check_open_unit('alpha', opts.alpha);
ivm_check_open_unit('op_tau', opts.op_tau);
ivm_check_open_unit('op_alpha', opts.op_alpha);
alpha = double(opts.alpha);
op_tau = double(opts.op_tau);
op_alpha = double(opts.op_alpha);
ar_ci = ivm_check_bool('ar_ci', opts.ar_ci);
lasso = ivm_check_bool('lasso', opts.lasso);
show_selection = ivm_check_bool('show_selection', opts.show_selection);
lasso_verbose = ivm_check_bool('lasso_verbose', opts.lasso_verbose);
lasso_penalty = ivm_check_choice('lasso_penalty', opts.lasso_penalty, {'plugin', 'cv'});
hw = opts.ar_grid_halfwidth;
gp = opts.ar_grid_points;
if ~(isnumeric(hw) && isscalar(hw) && isreal(hw)) || ~(isnumeric(gp) && isscalar(gp) && isreal(gp))
    error('ar_grid_halfwidth and ar_grid_points must be numbers.');
end
if ~(isfinite(hw) && hw > 0)
    error('ar_grid_halfwidth must be positive; got %s.', ivm_repr(hw));
end
if ~(gp == round(gp) && gp >= 3)
    error('ar_grid_points must be an integer >= 3; got %s.', ivm_repr(gp));
end
hw = double(hw);
gp = double(gp);
if lasso
    ivm_check_posnum('lasso_c_slack', opts.lasso_c_slack);
    if ~isempty(opts.lasso_alpha_p)
        ivm_check_open_unit('lasso_alpha_p', opts.lasso_alpha_p);
    end
    v = opts.lasso_n_iter_max;
    if ~(isnumeric(v) && isscalar(v) && isreal(v) && v == round(v) && v >= 1)
        error('lasso_n_iter_max must be a positive integer; got %s.', ivm_repr(v));
    end
    v = opts.lasso_tol;
    if ~(isnumeric(v) && isscalar(v) && isreal(v) && v > 0)
        error('lasso_tol must be positive; got %s.', ivm_repr(v));
    end
    v = opts.lasso_cv_folds;
    if ~(isnumeric(v) && isscalar(v) && isreal(v) && v == round(v) && v >= 2)
        error('lasso_cv_folds must be an integer >= 2; got %s.', ivm_repr(v));
    end
    v = opts.lasso_lambda_override;
    if ~isempty(v) && ~(isnumeric(v) && isscalar(v) && isreal(v) && v >= 0)
        error('lasso_lambda_override must be >= 0; got %s.', ivm_repr(v));
    end
end
expand_x = ivm_norm_text_list(opts.expand_x, 'expand_x');
expand_dummies = ivm_norm_text_list(opts.expand_dummies, 'expand_dummies');

% Spatial arguments are used only with se='spatial' (2026-10: an error otherwise).
clusterid = ivm_opt_name(opts.clusterid, 'clusterid');
weights = ivm_opt_name(opts.weights, 'weights');
lat = ivm_opt_name(opts.lat, 'lat');
lon = ivm_opt_name(opts.lon, 'lon');
spatial_cutoff = opts.spatial_cutoff;
spatial_kernel = opts.spatial_kernel;
spatial_metric = opts.spatial_metric;
given = {};
if ~isempty(lat), given{end+1} = 'lat'; end
if ~isempty(lon), given{end+1} = 'lon'; end
if ~isempty(spatial_cutoff), given{end+1} = 'spatial_cutoff'; end
if ~isempty(spatial_kernel), given{end+1} = 'spatial_kernel'; end
if ~isempty(spatial_metric), given{end+1} = 'spatial_metric'; end
if strcmp(se, 'spatial')
    if isempty(lat) || isempty(lon)
        error('se=''spatial'' requires lat and lon column names.');
    end
    if isempty(spatial_cutoff)
        error(['se=''spatial'' requires spatial_cutoff (distance bandwidth). ' ...
               'Conley (1999) results are sensitive to this choice; report ' ...
               'results at multiple cutoffs as a robustness check.']);
    end
    if ~(isnumeric(spatial_cutoff) && isscalar(spatial_cutoff) && isreal(spatial_cutoff))
        error('spatial_cutoff must be a positive number; got %s.', ivm_repr(spatial_cutoff));
    end
    if ~(isfinite(spatial_cutoff) && spatial_cutoff > 0)
        error('spatial_cutoff must be positive; got %s.', ivm_pystr(spatial_cutoff));
    end
    spatial_cutoff = double(spatial_cutoff);
    if isempty(spatial_kernel), spatial_kernel = 'uniform'; end
    if isempty(spatial_metric), spatial_metric = 'haversine'; end
    spatial_kernel = ivm_check_choice('spatial_kernel', spatial_kernel, {'uniform', 'bartlett'});
    spatial_metric = ivm_check_choice('spatial_metric', spatial_metric, {'haversine', 'euclidean'});
elseif ~isempty(given)
    error(['%s given but se=''%s'': the spatial arguments are used only ' ...
           'with se=''spatial''.'], strjoin(given, ', '), se);
end

% Cluster ids are used only with cluster SEs.
if strcmp(se, 'cluster') && isempty(clusterid)
    error('se=''cluster'' requires clusterid (the cluster variable).');
end
if ~isempty(clusterid) && ~strcmp(se, 'cluster')
    ivm_note(quiet, '  Note: clusterid ignored (se != ''cluster'').');
end

% -- Weight-type / SE-regime consistency --
if ~isempty(weights) && strcmp(wtype, 'pweight') && strcmp(se, 'classical')
    ivm_note(quiet, ['  Note: pweight + classical SE is inconsistent. ' ...
                     'Upgrading SE to ''robust'' (HC1 sandwich).']);
    se = 'robust';
elseif ~isempty(weights) && strcmp(wtype, 'aweight') && strcmp(se, 'classical')
    ivm_note(quiet, ['  Note: aweight + classical SE assumes the aweight ' ...
                     'specification is correctly calibrated. Consider ' ...
                     '''robust'' or ''cluster'' as a defensive default.']);
end
if point_only
    se_v = 'none';                 % variance regime actually computed
    ar_ci = false;
else
    se_v = se;
end

% -- Roles, missing values, listwise deletion, zero weights --
prep = ivm_prepare_sample(data, Y, T, M, Z, opts.X, opts.dummies, opts.fe, ...
                          se, clusterid, weights, lat, lon);
Y = prep.Y; T = prep.T; M = prep.M;
Z = prep.Z; X = prep.X; dummies = prep.dummies;
fe_cols = prep.fe_cols; clusterid = prep.clusterid; df = prep.df;
note_X = ivm_summarize_resolved('X', X);
note_D = ivm_summarize_resolved('dummies', dummies);
if ~isempty(note_X), ivm_note(quiet, note_X); end
if ~isempty(note_D), ivm_note(quiet, note_D); end
for i = 1:numel(prep.notes)
    ivm_note(quiet, prep.notes{i});
end
is_weighted = ~isempty(weights);
if prep.n_zero_weights > 0
    ivm_note(quiet, sprintf('  Note: dropping %d observation(s) with zero weight.', ...
                            prep.n_zero_weights));
end

% Expand dummy variables into indicator columns (on the estimation sample)
dummy_col_names = {};
if ~isempty(dummies)
    [df, dummy_col_names] = ivm_expand_dummies(df, dummies);
    if ~lasso
        X = [X, dummy_col_names];
    end
end

% Drop singleton FE groups (alternating passes until a full pass drops
% nothing; reghdfe convention). Group counts after convergence.
n_dropped_singletons = 0;
n_fe_groups = 0;
if ~isempty(fe_cols)
    while true
        dropped_this_pass = 0;
        for d = 1:numel(fe_cols)
            [df, ~, n_drop] = ivm_drop_singletons(df, fe_cols{d});
            dropped_this_pass = dropped_this_pass + n_drop;
        end
        n_dropped_singletons = n_dropped_singletons + dropped_this_pass;
        if dropped_this_pass == 0
            break;
        end
    end
    if n_dropped_singletons > 0
        ivm_note(quiet, sprintf(['  Note: dropped %d singleton observations ' ...
            'across %d FE dimension(s).'], n_dropped_singletons, numel(fe_cols)));
    end
end

n = ivm_df_nrows(df);
if n == 0
    error('No observations remaining after dropping missing values and singletons.');
end
n_clusters = [];
if ~isempty(clusterid)
    [~, n_clusters] = ivm_factorize(ivm_df_get(df, clusterid));
    if n_clusters < 2
        error('cluster SEs need at least 2 clusters.');
    end
end

% Absorbed FE degrees of freedom (rank of the FE dummies beyond the
% intercept); with cluster SEs, FE nested in the clusters do not count
% toward CR1 (reghdfe/fixest).
df_fe_absorbed = 0;
df_vcov = 0;
fe_nested_in_cluster = [];
if ~isempty(fe_cols)
    if strcmp(se, 'cluster')
        [df_fe_absorbed, df_fe_cr1, n_fe_groups, nested] = ivm_fe_df(df, fe_cols, clusterid);
        fe_nested_in_cluster = nested;
        df_vcov = df_fe_cr1;
    else
        [df_fe_absorbed, ~, n_fe_groups] = ivm_fe_df(df, fe_cols, []);
        df_vcov = df_fe_absorbed;
    end
end

% -- LASSO: build the rich expansion BEFORE FE demeaning --
rich_col_names = {};
rich_meta = [];
rich_in_df = {};
keep_idx_lasso = [];
if lasso
    [W_rich_levels, rich_col_names, rich_meta] = ivm_build_rich_expansion( ...
        df, X, dummy_col_names, expand_x, expand_dummies, opts.expand_override, quiet);
    rich_in_df = cell(1, numel(rich_col_names));
    for j = 1:numel(rich_col_names)
        nm = rich_col_names{j};
        if ivm_df_has(df, nm)
            rich_in_df{j} = sprintf('_RX_%d_%s', j - 1, nm);   % 0-based j, like Python
        else
            rich_in_df{j} = nm;
        end
    end
    for j = 1:numel(rich_in_df)
        df = ivm_df_set(df, rich_in_df{j}, W_rich_levels(:, j));
    end
    lasso_keep_resolved = ivm_resolve_lasso_keep(opts.lasso_keep, rich_col_names);
    keep_idx_lasso = zeros(1, numel(lasso_keep_resolved));
    for i = 1:numel(lasso_keep_resolved)
        keep_idx_lasso(i) = find(strcmp(rich_col_names, lasso_keep_resolved{i}), 1);
    end
    if ~isempty(rich_meta) && rich_meta.p > 0
        ivm_note(quiet, sprintf(['  Rich expansion built: p = %d ' ...
            '(originals_x=%d, originals_dummies=%d, sq=%d, cb=%d, ' ...
            'x%sx=%d, d%sx=%d, d%sd=%d)'], rich_meta.p, ...
            rich_meta.n_originals_x, rich_meta.n_originals_dummies, ...
            rich_meta.n_squares, rich_meta.n_cubes, ...
            ivm_uch(215), rich_meta.n_interactions_xx, ...
            ivm_uch(215), rich_meta.n_interactions_dx, ...
            ivm_uch(215), rich_meta.n_interactions_dd));
    end
    if ~isempty(lasso_keep_resolved)
        ivm_note(quiet, sprintf('  LASSO keep (%d): %s%s', ...
            numel(lasso_keep_resolved), ...
            strjoin(lasso_keep_resolved(1:min(8, numel(lasso_keep_resolved))), ', '), ...
            ivm_ellipsis(numel(lasso_keep_resolved) > 8)));
    end
    control_cols = rich_in_df;
else
    control_cols = X;
end

% Weights (all positive after the zero-weight drop)
if is_weighted
    w_vec = ivm_col_num(ivm_df_get(df, weights), weights);
else
    w_vec = ones(n, 1);
end
sqrt_w = sqrt(w_vec);

% Reference scales for the collinearity and identification checks: the
% weighted variation of each column around its weighted mean BEFORE the FE
% are absorbed.
ref_norm = @(v) sqrt(sum(w_vec .* (v - sum(w_vec .* v) / sum(w_vec)) .^ 2));
ref_ctrl = containers.Map('KeyType', 'char', 'ValueType', 'double');
for j = 1:numel(control_cols)
    ref_ctrl(control_cols{j}) = ref_norm(ivm_col_num(ivm_df_get(df, control_cols{j}), control_cols{j}));
end
role_names = [{T, M}, Z];
ref_role = containers.Map('KeyType', 'char', 'ValueType', 'double');
for j = 1:numel(role_names)
    ref_role(role_names{j}) = ref_norm(ivm_col_num(ivm_df_get(df, role_names{j}), role_names{j}));
end

% Demean within FE groups (weighted group means under weights)
analysis_cols = ivm_unique_stable([{Y, T, M}, Z, control_cols]);
if ~isempty(fe_cols)
    if numel(fe_cols) == 1
        df = ivm_demean_within_groups(df, analysis_cols, fe_cols{1}, w_vec);
    else
        df = ivm_demean_twoway(df, analysis_cols, fe_cols, w_vec, quiet);
    end
end

y = ivm_col_num(ivm_df_get(df, Y), Y);
t = ivm_col_num(ivm_df_get(df, T), T);
m = ivm_col_num(ivm_df_get(df, M), M);
z_mat = ivm_extract_mat(df, Z);
if ~isempty(control_cols)
    x_mat = ivm_extract_mat(df, control_cols);
else
    x_mat = [];
end
if ~isempty(clusterid)
    cluster_ids = ivm_df_get(df, clusterid);
else
    cluster_ids = [];
end

if strcmp(se, 'spatial')
    lat_arr = ivm_col_num(ivm_df_get(df, lat), lat);
    lon_arr = ivm_col_num(ivm_df_get(df, lon), lon);
    if strcmp(spatial_metric, 'haversine')
        if any(abs(lat_arr) > 90.0) || any(abs(lon_arr) > 180.0)
            error(['Lat/lon values are outside the (-90, 90) x (-180, 180) ' ...
                   'degree range. Check that you have lat/lon in degrees, or ' ...
                   'switch to spatial_metric=''euclidean'' for projected ' ...
                   'coordinates.']);
        end
    end
else
    lat_arr = [];
    lon_arr = [];
end

% -- SE plumbing shared by every 2SLS call in this analysis --
sekw = struct();
sekw.se = se_v;
sekw.cluster_gidx = [];
sekw.cluster_G = 0;
sekw.spatial_A = [];
spatial_A = [];
if strcmp(se, 'spatial') && (~point_only || lasso)
    [spatial_A, ~, nnz_count] = ivm_spatial_adjacency(lat_arr, lon_arr, ...
        spatial_cutoff, spatial_kernel, spatial_metric);
    share = nnz_count / n / n;
    if strcmp(spatial_kernel, 'uniform') && nnz_count >= n * n
        error(['the spatial cutoff covers all pairs; Conley SEs are degenerate ' ...
               '(0 by construction). Choose a smaller spatial_cutoff.']);
    end
    if share > 0.5
        ivm_note(quiet, sprintf(['  Warning: the spatial cutoff puts %.0f%% of all ' ...
            'pairs within the cutoff (average neighbor share > 50%%); Conley SEs ' ...
            'are unreliable in this regime. Consider a smaller cutoff.'], 100 * share));
    end
end
if strcmp(se_v, 'cluster')
    [sekw.cluster_gidx, sekw.cluster_G] = ivm_factorize(cluster_ids);
elseif strcmp(se_v, 'spatial')
    sekw.spatial_A = spatial_A;
end

y_raw = y; t_raw = t; m_raw = m; z_mat_raw = z_mat;
x_mat_raw = x_mat;

% -- LASSO step (raw-scale, FE-demeaned data) --
lasso_info = [];
if lasso
    if isempty(x_mat_raw) || size(x_mat_raw, 2) == 0
        lasso_info = struct();
        lasso_info.S = false(0, 1);
        lasso_info.S_Y = false(0, 1);
        lasso_info.S_T = false(0, 1);
        lasso_info.S_M = false(0, 1);
        lasso_info.S_Z = false(0, 1);
        lasso_info.S_Z_per_inst = {};
        details0 = struct();
        details0.Y = [];
        details0.T = [];
        details0.M = [];
        details0.Z = {};
        lasso_info.details = details0;
        lasso_info.union_names = {};
        lasso_info.n_selected_union = 0;
        lasso_info.n_total = 0;
        x_mat_raw = [];
        X = {};
        control_cols = {};
    else
        if ~isempty(cluster_ids)
            [cl_gidx, cl_G] = ivm_factorize(cluster_ids);
        else
            cl_gidx = [];
            cl_G = 0;
        end
        lopts = struct();
        lopts.n_iter_max = double(opts.lasso_n_iter_max);
        lopts.tol = double(opts.lasso_tol);
        lopts.c_slack = double(opts.lasso_c_slack);
        lopts.alpha_p = opts.lasso_alpha_p;
        lopts.lambda_method = lasso_penalty;
        lopts.lambda_override = opts.lasso_lambda_override;
        if strcmp(se, 'spatial')
            lopts.spatial_A = spatial_A;
        else
            lopts.spatial_A = [];
        end
        lopts.verbose = lasso_verbose;
        lasso_info = ivm_quadruple_lasso(y_raw, t_raw, m_raw, z_mat_raw, ...
            x_mat_raw, rich_col_names, w_vec, cl_gidx, cl_G, keep_idx_lasso, lopts);
        sel_idx = find(lasso_info.S(:).');
        if isempty(sel_idx)
            x_mat_raw = [];
            X = {};
            control_cols = {};
        else
            x_mat_raw = x_mat_raw(:, sel_idx);
            X = rich_col_names(sel_idx);
            control_cols = rich_in_df(sel_idx);
        end
        ivm_note(quiet, sprintf(['  LASSO selected %d of %d candidate ' ...
            'variables (%s penalty).'], lasso_info.n_selected_union, ...
            lasso_info.n_total, lasso_penalty));
    end
end

% -- Collinear controls: drop exactly redundant columns of [1, X] --
% Tested in order on the weighted, FE-demeaned design; a column whose
% residual on [1, earlier columns] is below 1e-9 of its variation before the
% FE were absorbed is dropped with a note.
dropped_collinear = {};
if ~isempty(x_mat_raw) && size(x_mat_raw, 2) > 0
    A = x_mat_raw .* sqrt_w;
    refs = zeros(1, numel(control_cols));
    for j = 1:numel(control_cols)
        refs(j) = ref_ctrl(control_cols{j});
    end
    [kept, dropped_collinear] = ivm_greedy_independent(A, X, 1e-9, sqrt_w, refs);
    if ~isempty(dropped_collinear)
        ivm_note(quiet, sprintf('  Note: dropped (collinear): %s', strjoin(dropped_collinear, ', ')));
        if isempty(kept)
            x_mat_raw = [];
        else
            x_mat_raw = x_mat_raw(:, kept);
        end
        X = X(kept);
        control_cols = control_cols(kept);
    end
end

% -- sqrt(w) premultiplication of every analysis array --
y = y_raw .* sqrt_w;
t = t_raw .* sqrt_w;
m = m_raw .* sqrt_w;
z_mat = z_mat_raw .* sqrt_w;
ones_v = sqrt_w;
if ~isempty(x_mat_raw)
    x_mat = x_mat_raw .* sqrt_w;
    exog = [ones_v, x_mat];
    Z_base = [ones_v, z_mat, x_mat];
else
    x_mat = [];
    exog = ones_v;
    Z_base = [ones_v, z_mat];
end

z_crit = ivm_norminv(1 - alpha / 2);
K_inst = size(z_mat, 2);
n_exog = size(exog, 2);

% -- Identification checks (paper Theorem 1 rank conditions, eq. 26) --
for j = 1:K_inst
    if ivm_rel_resid_var(z_mat(:, j), exog, ref_role(Z{j})) < 1e-12
        error(['instrument %s has no variation after the controls/fixed ' ...
               'effects (%s_ZZ = 0, paper eq. 26).'], Z{j}, ivm_uch(963));
    end
end
if K_inst > 1
    Zt_ = ivm_residualize(z_mat, exog);
    refs_z = zeros(1, K_inst);
    for j = 1:K_inst
        refs_z(j) = ref_role(Z{j});
    end
    [~, dz] = ivm_greedy_independent(Zt_, Z, 1e-9, [], refs_z);
    if ~isempty(dz)
        error('instruments are collinear after the controls/fixed effects (%s).', ...
              strjoin(dz, ', '));
    end
end
if ivm_rel_resid_var(t, exog, ref_role(T)) < 1e-12
    error('treatment %s has no variation after the controls/fixed effects.', T);
end
if ivm_rel_resid_var(m, exog, ref_role(M)) < 1e-12
    error('mediator %s has no variation after the controls/fixed effects.', M);
end
r_mt = ivm_residualize(m, [exog, t]);
if (r_mt' * r_mt) / ref_role(M)^2 < 1e-12
    error(['the mediator is collinear with the treatment after the ' ...
           'controls/fixed effects.']);
end

% Not enough observations for the largest regression ([1, X, Z, T])?
k_max = n_exog + K_inst + 1;
if n - k_max - df_fe_absorbed < 1
    error(['not enough observations: n = %d, but the largest regression has %d ' ...
           'coefficients plus %d absorbed fixed-effect degrees of freedom ' ...
           '(n must exceed their sum).'], n, k_max, df_fe_absorbed);
end

% ======================================================================
% 1. FIRST STAGE: T ~ Z + X (pi1, paper eq. 27; joint Wald F on the
%    excluded instruments under the same SE regime; OP effective F and the
%    simplified-test critical value at the estimated K_eff)
% ======================================================================
X_fs = Z_base;
res_fs = ivm_twosls_with_vcov(t, X_fs, [], sekw, df_vcov);
z_idx_fs = 2:(1 + K_inst);
beta_Z_fs = res_fs.coef(z_idx_fs);
V_ZZ_fs = res_fs.vcov(z_idx_fs, z_idx_fs);
pi1 = beta_Z_fs(:).';
se_pi1 = res_fs.se(z_idx_fs).';
F_TZ = ivm_wald_over_K(beta_Z_fs, V_ZZ_fs);
z_resid_fs = ivm_residualize(z_mat, exog);
Q_zz_fs = z_resid_fs' * z_resid_fs / n;
if K_inst == 1
    F_TZ_op = F_TZ;
elseif ivm_is_pd(V_ZZ_fs)
    F_TZ_op = ivm_op_effective_F(beta_Z_fs, V_ZZ_fs, Q_zz_fs);
else
    F_TZ_op = NaN;
end
[K_eff_TZ, F_TZ_op_critval] = ivm_op_keff_critval(V_ZZ_fs, Q_zz_fs, op_tau, op_alpha);

% ======================================================================
% 2. TREATMENT-MEDIATOR: M ~ T | Z + X (2SLS, instrument T with Z)
% ======================================================================
X_gamma = [exog, t];
res_gamma = ivm_twosls_with_vcov(m, X_gamma, Z_base, sekw, df_vcov);
idx_T_in_gamma = size(X_gamma, 2);
gamma_hat = res_gamma.coef(idx_T_in_gamma);
se_gamma = res_gamma.se(idx_T_in_gamma);

% ======================================================================
% 3. CONDITIONAL 2SLS: Y ~ M + T | Z + T + X (paper eqs. 32-34)
% ======================================================================
X_theta = [exog, m, t];
Z_theta = [Z_base, t];
res_theta = ivm_twosls_with_vcov(y, X_theta, Z_theta, sekw, df_vcov);
idx_M = n_exog + 1;
idx_T = n_exog + 2;
theta_hat = res_theta.coef(idx_M);
tau_hat = res_theta.coef(idx_T);
se_theta = res_theta.se(idx_M);
se_tau = res_theta.se(idx_T);

% ======================================================================
% 4. TOTAL EFFECT: Y ~ T | Z + X (standard 2SLS)
% ======================================================================
X_total = [exog, t];
res_total = ivm_twosls_with_vcov(y, X_total, Z_base, sekw, df_vcov);
idx_T_total = size(X_total, 2);
total_hat = res_total.coef(idx_T_total);
se_total = res_total.se(idx_T_total);

% ======================================================================
% 5. CONDITIONAL FIRST-STAGE F: M ~ Z + T + X (paper eq. 32, Remark 3)
% ======================================================================
X_cfs = [exog, z_mat, t];
res_cfs = ivm_twosls_with_vcov(m, X_cfs, [], sekw, df_vcov);
z_idx_cfs = (n_exog + 1):(n_exog + K_inst);
beta_Z_cfs = res_cfs.coef(z_idx_cfs);
V_ZZ_cfs = res_cfs.vcov(z_idx_cfs, z_idx_cfs);
F_MZ_T = ivm_wald_over_K(beta_Z_cfs, V_ZZ_cfs);
F_MZ_T_individual = NaN(1, K_inst);
for i = 1:K_inst
    ii = z_idx_cfs(i);
    if res_cfs.vcov(ii, ii) > 0
        F_MZ_T_individual(i) = res_cfs.coef(ii)^2 / res_cfs.vcov(ii, ii);
    end
end
z_resid_cfs = ivm_residualize(z_mat, [exog, t]);
Q_zz_cfs = z_resid_cfs' * z_resid_cfs / n;
if K_inst == 1
    F_MZ_T_op = F_MZ_T;
elseif ivm_is_pd(V_ZZ_cfs)
    F_MZ_T_op = ivm_op_effective_F(beta_Z_cfs, V_ZZ_cfs, Q_zz_cfs);
else
    F_MZ_T_op = NaN;
end
[K_eff_MZ_T, F_MZ_T_op_critval] = ivm_op_keff_critval(V_ZZ_cfs, Q_zz_cfs, op_tau, op_alpha);

% ======================================================================
% 6. DECOMPOSITION CHECK + delta-method SE for the indirect effect
% ======================================================================
indirect = theta_hat * gamma_hat;
decomp_check = tau_hat + indirect;
decomp_error = abs(decomp_check - total_hat);
decomp_error_rel = decomp_error / max([abs(total_hat), abs(tau_hat) + abs(indirect), 1e-300]);

cov_gamma_theta = ivm_cross_equation_cov(res_gamma, res_theta, idx_T_in_gamma, ...
                                         idx_M, sekw, df_vcov);
var_indirect = gamma_hat^2 * se_theta^2 + theta_hat^2 * se_gamma^2 ...
               + 2.0 * theta_hat * gamma_hat * cov_gamma_theta;
if var_indirect < 0
    if var_indirect > -1e-8 * max((gamma_hat * se_theta)^2, 1e-300)
        var_indirect = 0.0;
    else
        ivm_note(quiet, sprintf(['  Warning: delta-method variance for the ' ...
            'indirect effect is negative (%.2e). This typically reflects a ' ...
            'poorly identified first stage. Reporting NaN; consider ' ...
            'ivmediate_bootstrap for inference.'], var_indirect));
        var_indirect = NaN;
    end
end
se_indirect = ivm_sqrtnan(var_indirect);

% ======================================================================
% 7. KAPPA-BASED SENSITIVITY (paper Section 5.1.1, Proposition 4; Hajek-
%    form weighted moments on the raw arrays, X partialed out by WLS)
% ======================================================================
if ~isempty(x_mat_raw) && size(x_mat_raw, 2) > 0
    X_partial = [ones(n, 1), x_mat_raw];
else
    X_partial = ones(n, 1);
end
fit_t = ivm_wls_fit(X_partial, t_raw, w_vec);
t_orth = fit_t.residuals;
fit_m = ivm_wls_fit(X_partial, m_raw, w_vec);
m_orth = fit_m.residuals;
z_mat_orth = zeros(n, K_inst);
for k = 1:K_inst
    fit_z = ivm_wls_fit(X_partial, z_mat_raw(:, k), w_vec);
    z_mat_orth(:, k) = fit_z.residuals;
end
fs_kappa = ivm_wls_fit([ones(n, 1), z_mat_orth], t_orth, w_vec);
eps_T_hat = fs_kappa.residuals;
eps_M_hat = m_orth - gamma_hat * t_orth;
eps_Y_hat = res_theta.resid ./ sqrt_w;

c_TM_hat = ivm_wcov(eps_T_hat, eps_M_hat, w_vec);
sigma_eT_hat = ivm_wstd(eps_T_hat, w_vec);
sigma_eY_hat = ivm_wstd(eps_Y_hat, w_vec);

kappa_star = theta_hat * c_TM_hat;
% Inference breakeven (D1): kappa*_CI = (theta - sign(theta) z SE) c_TM;
% uninformative (flagged) when the CI already contains 0.
kappa_star_CI = (theta_hat - sign(theta_hat) * z_crit * se_theta) * c_TM_hat;
if isfinite(se_theta)
    kappa_star_CI_informative = abs(theta_hat) > z_crit * se_theta;
else
    kappa_star_CI_informative = [];
end
if sigma_eT_hat > 0 && sigma_eY_hat > 0
    rho_TY_star = kappa_star / (sigma_eT_hat * sigma_eY_hat);
else
    rho_TY_star = NaN;
end

kappa_frac_grid = [0.0, 0.10, 0.25, 0.50, 0.75, 0.90, 1.00];
if isnan(kappa_star) || kappa_star == 0
    kappa_grid = 0.0;
    kappa_frac_used = 0.0;
else
    kappa_grid = kappa_frac_grid * kappa_star;
    kappa_frac_used = kappa_frac_grid;
end
kappa_adj_table = struct('kappa', {}, 'kappa_frac', {}, 'theta_adj', {}, ...
                         'ci_lo', {}, 'ci_hi', {});
for i = 1:numel(kappa_grid)
    th_adj = theta_hat - kappa_grid(i) / c_TM_hat;
    row = struct('kappa', kappa_grid(i), 'kappa_frac', kappa_frac_used(i), ...
                 'theta_adj', th_adj, 'ci_lo', th_adj - z_crit * se_theta, ...
                 'ci_hi', th_adj + z_crit * se_theta);
    kappa_adj_table(end+1) = row; %#ok<AGROW>
end
% Strength of the conditional first stage (App. F.4 step 3): true iff the
% effective F_MZ|T reaches the OP critical value at the user's op_tau; []
% when either is not finite. ivmediate_print notes it below the kappa table.
if isfinite(F_MZ_T_op) && isfinite(F_MZ_T_op_critval)
    F_MZ_T_above_crit = F_MZ_T_op >= F_MZ_T_op_critval;
else
    F_MZ_T_above_crit = [];
end

% ======================================================================
% 8. R^2-BASED SENSITIVITY (paper App. F.2-F.3, eqs. 111-125)
% ======================================================================
W_fwl = [exog, t];
fit_y = ivm_lm_fit(W_fwl, y);
Y_tilde = fit_y.residuals;
fit_mt = ivm_lm_fit(W_fwl, m);
M_tilde = fit_mt.residuals;
fit_zt = ivm_lm_fit(W_fwl, z_mat(:, 1));
Z_tilde = fit_zt.residuals;

X_Z = [ones_v, Z_tilde];
df_adj_rv = (size(W_fwl, 2) - 1) + df_vcov;

res_rf = ivm_twosls_with_vcov(Y_tilde, X_Z, [], sekw, df_adj_rv);
t_rf = res_rf.coef(2) / res_rf.se(2);
F_reduced_form = t_rf^2;

res_fs_cond = ivm_twosls_with_vcov(M_tilde, X_Z, [], sekw, df_adj_rv);
t_fs_cond = res_fs_cond.coef(2) / res_fs_cond.se(2);
F_first_stage_cond = t_fs_cond^2;

% Residual df of the regression of Y (or M) on [1, X, T, Z] (sensemakr's
% model$df.residual).
df_rv = n - size(W_fwl, 2) - 1 - df_fe_absorbed;
t_star_rv = abs(ivm_tinv(alpha / 2.0, df_rv - 1));
f_star_rv = t_star_rv / sqrt(df_rv - 1);
if isfinite(F_reduced_form) && isfinite(F_first_stage_cond)
    xrv_rf = ivm_xrv_ols(F_reduced_form, df_rv, f_star_rv);
    xrv_fs = ivm_xrv_ols(F_first_stage_cond, df_rv, f_star_rv);
    rv_rf = ivm_rv_ols(F_reduced_form, df_rv, f_star_rv);
    rv_fs = ivm_rv_ols(F_first_stage_cond, df_rv, f_star_rv);
    XRV_theta = min(xrv_rf, xrv_fs);
    RV_theta = min(rv_rf, rv_fs);
else
    xrv_rf = NaN; xrv_fs = NaN; rv_rf = NaN; rv_fs = NaN;
    XRV_theta = NaN; RV_theta = NaN;
end

% ======================================================================
% 9. CONFIDENCE INTERVALS + ANDERSON-RUBIN SETS
% ======================================================================
ci_theta = [theta_hat - z_crit * se_theta, theta_hat + z_crit * se_theta];
ci_tau = [tau_hat - z_crit * se_tau, tau_hat + z_crit * se_tau];
ci_gamma = [gamma_hat - z_crit * se_gamma, gamma_hat + z_crit * se_gamma];
ci_total = [total_hat - z_crit * se_total, total_hat + z_crit * se_total];
if isfinite(se_indirect)
    ci_indirect = [indirect - z_crit * se_indirect, indirect + z_crit * se_indirect];
else
    ci_indirect = [NaN, NaN];
end

ar_theta = []; ar_gamma = []; ar_total = [];
if ar_ci
    W_ar_theta = [exog, t, z_mat];
    z_cols_theta = (n_exog + 2):(n_exog + 1 + K_inst);
    ar_theta = ivm_ar_confidence_set(y, m, W_ar_theta, z_cols_theta, sekw, ...
        df_vcov, alpha, theta_hat, se_theta, hw, gp);
    z_cols_zb = 2:(1 + K_inst);
    ar_gamma = ivm_ar_confidence_set(m, t, Z_base, z_cols_zb, sekw, ...
        df_vcov, alpha, gamma_hat, se_gamma, hw, gp);
    ar_total = ivm_ar_confidence_set(y, t, Z_base, z_cols_zb, sekw, ...
        df_vcov, alpha, total_hat, se_total, hw, gp);
end

% Negative spatial variances: one warning per analysis, and only for the
% reported quantities.
if strcmp(se, 'spatial') && ~point_only
    rep_var = [res_gamma.vcov(idx_T_in_gamma, idx_T_in_gamma), ...
               res_theta.vcov(idx_M, idx_M), res_theta.vcov(idx_T, idx_T), ...
               res_total.vcov(idx_T_total, idx_T_total), ...
               res_rf.vcov(2, 2), res_fs_cond.vcov(2, 2)];
    if any(rep_var < 0) || ~(ivm_is_pd(V_ZZ_fs) && ivm_is_pd(V_ZZ_cfs))
        ivm_note(quiet, ['  Warning: the spatial (Conley) variance was negative ' ...
            'for some reported coefficients (the kernel matrix is not positive ' ...
            'semi-definite); their SEs and the statistics that use them are ' ...
            'reported as NaN.']);
    end
end

% ======================================================================
% 10. ASSEMBLE RESULTS (None in Python = [] here)
% ======================================================================
results = struct();
results.version = ivm_version();
results.n = n;
results.n_instruments = K_inst;
results.n_covariates = numel(X);
if ~isempty(dummies)
    results.n_dummies = numel(dummy_col_names);
else
    results.n_dummies = 0;
end
results.se_type = se;
results.alpha = alpha;
if isempty(fe_cols)
    results.fe = [];
else
    results.fe = fe_cols;
end
results.n_fe_groups = n_fe_groups;
results.n_dropped_singletons = n_dropped_singletons;
results.df_fe_absorbed = df_fe_absorbed;
results.fe_nested_in_cluster = fe_nested_in_cluster;
results.n_clusters = n_clusters;
results.dropped_collinear = dropped_collinear;
results.weighted = is_weighted;
if is_weighted
    results.wtype = wtype;
else
    results.wtype = [];
end
if strcmp(se, 'spatial')
    results.spatial_cutoff = spatial_cutoff;
    results.spatial_kernel = spatial_kernel;
    results.spatial_metric = spatial_metric;
else
    results.spatial_cutoff = [];
    results.spatial_kernel = [];
    results.spatial_metric = [];
end
results.theta_hat = theta_hat; results.se_theta = se_theta; results.ci_theta = ci_theta;
results.tau_hat = tau_hat; results.se_tau = se_tau; results.ci_tau = ci_tau;
results.gamma_hat = gamma_hat; results.se_gamma = se_gamma; results.ci_gamma = ci_gamma;
results.total_hat = total_hat; results.se_total = se_total; results.ci_total = ci_total;
results.indirect = indirect; results.se_indirect = se_indirect; results.ci_indirect = ci_indirect;
results.cov_gamma_theta = cov_gamma_theta;
results.pi1 = pi1;                 % 1 x K (a scalar for K = 1)
results.se_pi1 = se_pi1;
results.decomp_check = decomp_check;
results.decomp_error = decomp_error;
results.decomp_error_rel = decomp_error_rel;
results.F_TZ = F_TZ;
results.F_MZ_T = F_MZ_T;
results.F_MZ_T_individual = F_MZ_T_individual;
results.F_TZ_op = F_TZ_op;
results.F_MZ_T_op = F_MZ_T_op;
results.K_eff_TZ = K_eff_TZ;
results.K_eff_MZ_T = K_eff_MZ_T;
results.F_TZ_op_critval = F_TZ_op_critval;
results.F_MZ_T_op_critval = F_MZ_T_op_critval;
results.op_tau = op_tau;
results.op_alpha = op_alpha;
results.ar_ci = ar_ci;
results.ar_theta = ar_theta;
results.ar_gamma = ar_gamma;
results.ar_total = ar_total;
results.c_TM_hat = c_TM_hat;
results.rho_TM_hat = c_TM_hat;     % deprecated alias (pre-2026-07 name)
results.sigma_eT = sigma_eT_hat;
results.sigma_eY = sigma_eY_hat;
results.kappa_star = kappa_star;
results.kappa_star_CI = kappa_star_CI;
results.kappa_star_CI_informative = kappa_star_CI_informative;
results.rho_TY_star = rho_TY_star;
results.kappa_adj_table = kappa_adj_table;
results.F_MZ_T_above_crit = F_MZ_T_above_crit;
results.F_reduced_form = F_reduced_form;
results.F_first_stage_cond = F_first_stage_cond;
results.f_reduced_form = F_reduced_form;          % deprecated alias (= t^2)
results.f_first_stage_cond = F_first_stage_cond;  % deprecated alias (= t^2)
results.RV_theta = RV_theta;
results.XRV_theta = XRV_theta;
results.rv_rf = rv_rf;
results.rv_fs = rv_fs;
results.xrv_rf = xrv_rf;
results.xrv_fs = xrv_fs;
results.lasso = lasso;
if lasso
    results.lasso_penalty = lasso_penalty;
else
    results.lasso_penalty = [];
end
if lasso && strcmp(lasso_penalty, 'cv')
    results.lasso_seed = opts.lasso_seed;
else
    results.lasso_seed = [];
end
results.lasso_info = lasso_info;
results.lasso_expansion_meta = rich_meta;
results.lasso_show_selection = lasso && show_selection;
if lasso
    results.rich_col_names = rich_col_names;
else
    results.rich_col_names = [];
end
results.eps_T_hat = eps_T_hat;
results.eps_M_hat = eps_M_hat;
results.eps_Y_hat = eps_Y_hat;
end

function r = ivm_rel_resid_var(v_w, exog, ref)
% Residual variance of v on exog relative to the squared reference scale.
if ~(ref > 0)
    r = 0.0;
    return;
end
rr = ivm_residualize(v_w, exog);
r = (rr' * rr) / ref^2;
end

function F = ivm_wald_over_K(b, V)
% NaN unless the variance block is positive definite.
if ~ivm_is_pd(V)
    F = NaN;
    return;
end
F = (b(:)' * (V \ b(:))) / numel(b);
end

function v = ivm_xrv_ols(f_stat, df_rv, f_star_rv)
% XRV_{1,alpha} for an individual OLS coefficient (f_stat = t^2).
f_beta2 = abs(f_stat) / df_rv;
v = max((f_beta2 - f_star_rv^2) / (1.0 + f_beta2), 0.0);
end

function v = ivm_rv_ols(f_stat, df_rv, f_star_rv)
% RV_{1,alpha} for an individual OLS coefficient (f_stat = t^2).
f_beta = sqrt(abs(f_stat) / df_rv);
fqa = f_beta - f_star_rv;
if fqa < 0.0
    v = 0.0;
    return;
end
if f_beta < 1.0 / f_star_rv
    v = 0.5 * (sqrt(fqa^4 + 4.0 * fqa^2) - fqa^2);
    return;
end
v = ivm_xrv_ols(f_stat, df_rv, f_star_rv);
end

% ============================================================================
%  Print functions (every line goes through ivm_out; option 'ascii')
% ============================================================================

function cuts = ivm_validate_star_cutoffs(cutoffs)
if ~isnumeric(cutoffs) || ~isreal(cutoffs)
    error('star_cutoffs must be an iterable of 3 numbers; got %s.', ivm_repr(cutoffs));
end
if numel(cutoffs) ~= 3
    error('star_cutoffs must have exactly 3 entries; got %d.', numel(cutoffs));
end
cuts = double(cutoffs(:).');
if ~all(isfinite(cuts) & cuts > 0 & cuts < 1)
    error('star_cutoffs must all be strictly between 0 and 1; got %s.', mat2str(cuts));
end
cuts = sort(cuts);
end

function p = ivm_wald_pvalue(est, se)
if isempty(est) || isempty(se) || ~isfinite(est) || ~isfinite(se) || se <= 0
    p = NaN;
    return;
end
p = 2.0 * ivm_normsf(abs(est / se));
end

function s = ivm_sig_stars(pval, cutoffs)
if isempty(pval) || ~isfinite(pval)
    s = '';
elseif pval < cutoffs(1)
    s = '***';
elseif pval < cutoffs(2)
    s = '**';
elseif pval < cutoffs(3)
    s = '*';
else
    s = '';
end
end

function s = ivm_format_cut(c)
s = sprintf('%g', c);
end

function s = ivm_pyrepr(v)
% Python repr() of a float: shortest round-trip form, '.0' for integers.
if isempty(v)
    s = 'None';
elseif ~isfinite(v)
    if isnan(v)
        s = 'nan';
    elseif v > 0
        s = 'inf';
    else
        s = '-inf';
    end
elseif v == fix(v) && abs(v) < 1e16
    s = sprintf('%.1f', v);
else
    for p = 1:17
        s = sprintf('%.*g', p, v);
        if str2double(s) == v
            break;
        end
    end
    % Python switches to exponent form outside [1e-4, 1e16)
    if abs(v) >= 1e-4 && any(s == 'e')
        s = sprintf('%.*f', max(0, p - 1 - floor(log10(abs(v)))), v);
    end
end
end

function ivm_print_results(r, varargin)
%IVM_PRINT_RESULTS Formatted summary (mirrors the Python print_results).
% Options: 'title' (text), 'boot' (ivmediate_bootstrap output or []),
%          'stars' (logical), 'star_cutoffs' (3-vector), 'kappa_table'
%          (false; true also prints the kappa-adjusted estimates),
%          'r2_sensitivity' (false; true also prints the R^2-based
%          sensitivity section), 'ascii' ([] = auto).
defaults = struct();
defaults.title = 'IV Mediation Analysis';
defaults.boot = [];
defaults.stars = true;
defaults.star_cutoffs = [0.01, 0.05, 0.10];
defaults.kappa_table = false;
defaults.r2_sensitivity = false;
defaults.ascii = [];
popts = ivm_parse_opts(defaults, varargin, 'ivmediate_print');
[old_ascii, wstate] = ivm_enter(popts.ascii);
cleanup = onCleanup(@() ivm_leave(old_ascii, wstate)); %#ok<NASGU>
ivm_print_results_impl(r, popts);
end

function ivm_print_results_impl(r, popts)
title_str = ivm_norm_text_scalar(popts.title, 'title');
boot = popts.boot;
stars = ivm_check_bool('stars', popts.stars);
kappa_table = ivm_check_bool('kappa_table', popts.kappa_table);
r2_sensitivity = ivm_check_bool('r2_sensitivity', popts.r2_sensitivity);
star_cutoffs = popts.star_cutoffs;
if stars
    star_cutoffs = ivm_validate_star_cutoffs(star_cutoffs);
end
lvl = ivm_level_label(r.alpha);
using_boot = ~isempty(boot);
z_crit = ivm_norminv(1 - r.alpha / 2);
asc = ivm_ascii_state();

TH = ivm_uch(952); GA = ivm_uch(947); TAU = ivm_uch(964); KAP = ivm_uch(954);
RHO = ivm_uch(961); EPS = ivm_uch(949); ARR = ivm_uch(8594);
HBAR = ivm_uch(9472); NDASH = ivm_uch(8211); MDASH = ivm_uch(8212);
CUP = ivm_uch(8746);
SUP2 = ivm_uch(178); LAM = ivm_uch(955); PI_ = ivm_uch(960);
if asc
    TG = 'theta*gamma';     % as in Stata
else
    TG = [TH, GA];
end

ivm_out('\n%s\n', repmat('=', 1, 70));
ivm_out('  %s\n', title_str);
ivm_out('%s\n', repmat('=', 1, 70));
info_line = sprintf('  n = %s  |  Instruments: %d  |  Covariates: %d  |  SE: %s', ...
    ivm_thousands(r.n), r.n_instruments, r.n_covariates, r.se_type);
if isfield(r, 'n_clusters') && ~isempty(r.n_clusters)
    info_line = [info_line, sprintf(' (G = %d clusters)', r.n_clusters)];
end
if using_boot
    info_line = [info_line, sprintf('  +  Bootstrap (B=%d, %s)', boot.B, boot.method)];
end
ivm_out('%s\n', info_line);
if r.weighted
    ivm_out('  Weighted: %s\n', r.wtype);
end
if r.n_dummies > 0
    ivm_out('  Dummy indicators: %d\n', r.n_dummies);
end
if isfield(r, 'dropped_collinear') && ~isempty(r.dropped_collinear)
    ivm_out('  Dropped (collinear): %s\n', strjoin(r.dropped_collinear, ', '));
end
if ~isempty(r.fe)
    line = sprintf('  FE: %s  |  Groups: %d  |  Singletons dropped: %d  |  df absorbed: %d', ...
        strjoin(r.fe, ', '), r.n_fe_groups, r.n_dropped_singletons, r.df_fe_absorbed);
    nested = [];
    if isfield(r, 'fe_nested_in_cluster')
        nested = r.fe_nested_in_cluster;
    end
    if ~isempty(nested) && any(nested)
        line = [line, '  |  nested in clusters: ', strjoin(r.fe(logical(nested)), ', '), ...
                ' (not counted in the CR1 df)'];
    end
    ivm_out('%s\n', line);
end
if strcmp(r.se_type, 'spatial')
    ivm_out('  Spatial: cutoff = %s (metric=%s, kernel=%s)\n', ...
        ivm_pyrepr(r.spatial_cutoff), r.spatial_metric, r.spatial_kernel);
end
if r.lasso
    li = r.lasso_info;
    if isempty(li)
        n_sel = 0; n_tot = 0;
    else
        n_sel = li.n_selected_union; n_tot = li.n_total;
    end
    ivm_out('  LASSO: post-double selection, penalty=''%s'', selected %d/%d candidates\n', ...
        r.lasso_penalty, n_sel, n_tot);
end
ivm_out('%s\n', repmat(HBAR, 1, 70));

% -- Coefficient table (labels padded to 20 as in Python; ASCII mode: to the
%    longest ASCII label, the header and the rule widened to match) --
ivm_out('\n  MEDIATION ESTIMATES\n');
ci_header = sprintf('[%s CI]', lvl);
if using_boot
    se_header = 'Boot SE';
    ci_header = sprintf('[Boot %s CI]', lvl);
else
    se_header = 'SE';
end
row_names = {sprintf('%s_total (Wald)', TAU), ...
             sprintf('%s (direct T%sY)', TAU, ARR), ...
             sprintf('%s (mediator%sY)', TH, ARR), ...
             sprintf('%s (T%smediator)', GA, ARR), ...
             sprintf('%s (indirect)', TG)};
wl = ivm_label_width(row_names, 20);
ivm_out('  %s %s %s %s\n', ivm_pad_right('Parameter', wl), ...
    ivm_pad_left('Estimate', 10), ivm_pad_left(se_header, 10), ...
    ivm_pad_left(ci_header, 24));
ivm_out('  %s\n', repmat(HBAR, 1, 64 + wl - 20));
ests = [r.total_hat, r.tau_hat, r.theta_hat, r.gamma_hat, r.indirect];
if using_boot
    ses = [boot.boot_se.total_hat, boot.boot_se.tau_hat, boot.boot_se.theta_hat, ...
           boot.boot_se.gamma_hat, boot.boot_se.indirect_hat];
    cis = [boot.boot_ci.total_hat; boot.boot_ci.tau_hat; boot.boot_ci.theta_hat; ...
           boot.boot_ci.gamma_hat; boot.boot_ci.indirect_hat];
else
    ses = [r.se_total, r.se_tau, r.se_theta, r.se_gamma, r.se_indirect];
    cis = [r.ci_total; r.ci_tau; r.ci_theta; r.ci_gamma; r.ci_indirect];
end
for i = 1:5
    if isfinite(ses(i))
        if stars
            star_str = ivm_sig_stars(ivm_wald_pvalue(ests(i), ses(i)), star_cutoffs);
        else
            star_str = '';
        end
        if ~isempty(star_str)
            suffix = [' ', star_str];
        else
            suffix = '';
        end
        ivm_out('  %s %10.4f %10.4f   [%9.4f, %9.4f]%s\n', ...
            ivm_pad_right(row_names{i}, wl), ests(i), ses(i), cis(i, 1), cis(i, 2), suffix);
    else
        ivm_out('  %s %10.4f\n', ivm_pad_right(row_names{i}, wl), ests(i));
    end
end
if stars
    ivm_out('  Sig. codes: *** p<%s, ** p<%s, * p<%s\n', ...
        ivm_format_cut(star_cutoffs(1)), ivm_format_cut(star_cutoffs(2)), ...
        ivm_format_cut(star_cutoffs(3)));
end

% -- Decomposition check (silent unless violated): exact for K = 1 (relative
%    error above 1e-8 is an inconsistency); asymptotic for K > 1 (1e-2) --
if isfield(r, 'decomp_error_rel') && ~isempty(r.decomp_error_rel)
    rel = r.decomp_error_rel;
else
    rel = r.decomp_error / max([abs(r.total_hat), abs(r.tau_hat) + abs(r.indirect), 1e-300]);
end
if r.n_instruments == 1
    if rel > 1e-8
        ivm_out('\n  WARNING: decomposition identity violated.\n');
        ivm_out('  %s + %s = %.6f  vs  %s_total = %.6f\n', TAU, TG, ...
            r.decomp_check, TAU, r.total_hat);
        ivm_out('  relative error = %.2e  (threshold = 1e-08)\n', rel);
        ivm_out(['  For a single instrument this identity is exact; a nonzero ' ...
                 'error indicates an internal inconsistency. Please report.\n']);
    end
else
    if rel > 1e-2
        ivm_out('\n  Note: with K>1 instruments, %s + %s = %.4f and the Wald\n', ...
            TAU, TG, r.decomp_check);
        ivm_out('  %s_total = %.4f differ by %.1e (relative). This is expected %s the\n', ...
            TAU, r.total_hat, rel, MDASH);
        ivm_out('  objects are estimated equation-by-equation (not by joint GMM) and coincide\n');
        ivm_out('  only asymptotically. Not an error.\n');
    end
end

% -- First-stage diagnostics (the labels end at column 22, as the F lines) --
ivm_out('\n  FIRST-STAGE DIAGNOSTICS  (SE: %s)\n', r.se_type);
if isfield(r, 'pi1') && ~isempty(r.pi1)
    % K > 1: the label column is two characters wider (the widest label is
    % 'F_{M, Z|T, X} (OP eff.)')
    pi_lab = ivm_pad_right_dw(sprintf('%s_1 (T on Z | X)', PI_), 22 + 2 * (numel(r.pi1) > 1));
    if numel(r.pi1) > 1
        parts = cell(1, numel(r.pi1));
        for i = 1:numel(r.pi1)
            parts{i} = sprintf('%.4f (%.4f)', r.pi1(i), r.se_pi1(i));
        end
        ivm_out('  %s= %s\n', pi_lab, strjoin(parts, ', '));
    else
        ivm_out('  %s= %.4f (%.4f)\n', pi_lab, r.pi1, r.se_pi1);
    end
end
op_str_TZ = ivm_op_str(r.F_TZ_op_critval, ivm_getfield(r, 'K_eff_TZ', NaN), r.op_tau, r.op_alpha);
if r.n_instruments > 1
    ivm_out('  F_{T, Z|X}              = %.2f  (joint Wald)\n', r.F_TZ);
    ivm_out('  F_{T, Z|X} (OP eff.)    = %.2f%s\n', r.F_TZ_op, op_str_TZ);
else
    ivm_out('  F_{T, Z|X}            = %.2f%s\n', r.F_TZ, op_str_TZ);
end
op_str_MZ = ivm_op_str(r.F_MZ_T_op_critval, ivm_getfield(r, 'K_eff_MZ_T', NaN), r.op_tau, r.op_alpha);
if r.n_instruments > 1
    ivm_out('  F_{M, Z|T, X}           = %.2f  (joint Wald)\n', r.F_MZ_T);
    ivm_out('  F_{M, Z|T, X} (OP eff.) = %.2f%s\n', r.F_MZ_T_op, op_str_MZ);
    parts = cell(1, numel(r.F_MZ_T_individual));
    for i = 1:numel(r.F_MZ_T_individual)
        parts{i} = sprintf('%.1f', r.F_MZ_T_individual(i));
    end
    ivm_out('%s%s\n', ivm_pad_right_dw(sprintf('    per-instrument t%s:', SUP2), 28), strjoin(parts, '  '));
else
    ivm_out('  F_{M, Z|T, X}         = %.2f%s\n', r.F_MZ_T, op_str_MZ);
end

% -- Anderson-Rubin confidence sets (labels padded to 14 as in Python) --
ar_shown = r.ar_ci && ~isempty(r.ar_theta);
if ar_shown
    ivm_out('\n  ANDERSON%sRUBIN %s CONFIDENCE SETS  (weak-IV-robust, SE: %s)\n', ...
        NDASH, lvl, r.se_type);
    any_unbounded = false;
    labels = {sprintf('%s (M %s Y | T)', TH, ARR), ...
              sprintf('%s (T %s M)', GA, ARR), ...
              sprintf('%s_total', TAU)};
    wa = ivm_label_width(labels, 14);
    ars = {r.ar_theta, r.ar_gamma, r.ar_total};
    for i = 1:3
        ar = ars{i};
        if isempty(ar)
            continue;
        end
        ivm_out('  %s %s\n', ivm_pad_right(labels{i}, wa), ivm_format_ar_set(ar));
        if isfinite(ar.pval_at_zero)
            pv0_s = sprintf('%.4f', ar.pval_at_zero);
        else
            pv0_s = 'n/a';
        end
        if strcmp(ar.method, 'fieller')
            meth = 'exact Fieller inversion';
        else
            meth = 'exact statistic, grid search';
        end
        ivm_out('  %s   AR p-value at 0 = %s   (%s)\n', repmat(' ', 1, wa), pv0_s, meth);
        if ivm_ar_unbounded_or_disconnected(ar)
            any_unbounded = true;
        end
    end
    if any_unbounded
        ivm_out('  Note: an unbounded or disconnected AR set is the expected outcome under weak\n');
        ivm_out('  identification %s the data cannot bound the parameter at this confidence level.\n', MDASH);
    end
    ivm_out('  AR sets remain valid regardless of first-stage strength; compare with the Wald\n');
    ivm_out('  CIs above, which require the corresponding effective F to clear the OP threshold.\n');
end

% -- Kappa-based sensitivity: c_TM, kappa*, kappa*_CI and rho*_TY; the
%    kappa-adjusted table only on request (kappa_table) --
ivm_out('\n  %s-BASED SENSITIVITY (MCA diagnostics)\n', KAP);
if r.weighted
    ivm_out('  (Hajek-weighted moments, wtype = %s)\n', r.wtype);
end
if using_boot
    ivm_out('  (%s-frame CIs use bootstrap SE on %s)\n', KAP, TH);
end
if using_boot
    se_theta_use = boot.boot_se.theta_hat;
else
    se_theta_use = r.se_theta;
end
th = r.theta_hat;
kappa_star_CI_use = (th - sign(th) * z_crit * se_theta_use) * r.c_TM_hat;
informative = isfinite(se_theta_use) && abs(th) > z_crit * se_theta_use;
if informative
    flag = '';
else
    flag = '  (uninformative: the CI already contains 0)';
end
rho = r.rho_TY_star;
if isfinite(rho) && abs(rho) > 1
    rho_note = sprintf('  (|%s*| > 1: no feasible error correlation overturns %s)', ...
                       RHO, TH);
else
    rho_note = '';
end
klab = {sprintf('c_TM = Cov(%s_T, %s_M)', EPS, EPS), ...
        sprintf('%s* (breakeven)', KAP), ...
        sprintf('%s*_CI (inf. break.)', KAP), ...
        sprintf('%s*_TY (as corr.)', RHO)};
if asc
    wk = ivm_label_width(klab, 0);
    for i = 1:4
        klab{i} = [ivm_pad_right(klab{i}, wk), ' '];
    end
else
    gaps = {' ', '       ', '  ', '     '};          % Python's spacing
    for i = 1:4
        klab{i} = [klab{i}, gaps{i}];
    end
end
ivm_out('  %s= %.4f\n', klab{1}, r.c_TM_hat);
ivm_out('  %s= %.4f\n', klab{2}, r.kappa_star);
ivm_out('  %s= %.4f%s\n', klab{3}, kappa_star_CI_use, flag);
ivm_out('  %s= %.4f%s\n', klab{4}, rho, rho_note);
if using_boot
    keys = {'kappa_star', 'rho_TY_star'};
    labs = {[KAP, '*'], [RHO, '*_TY']};
    if asc
        labs = {ivm_pad_right(labs{1}, ivm_label_width(labs, 0)), ...
                ivm_pad_right(labs{2}, ivm_label_width(labs, 0))};
    end
    for i = 1:2
        if isfield(boot.boot_ci, keys{i})
            ci_b = boot.boot_ci.(keys{i});
            ivm_out('  %s bootstrap %s CI (%s) = [%.4f, %.4f]\n', labs{i}, ...
                ivm_level_label(boot.alpha), boot.ci_type, ci_b(1), ci_b(2));
        end
    end
end
if kappa_table
    ivm_out('\n  %s-adjusted estimates  (%s_0 swept as fraction of breakeven %s*):\n', ...
        KAP, KAP, KAP);
    hdr = {['% ', KAP, '*'], [KAP, '_0'], [TH, '_adj'], 'CI_lo', 'CI_hi'};
    wc = [6, 12, 10, 10, 10];
    if asc
        for j = 1:5
            wc(j) = max(wc(j), ivm_strwidth(hdr{j}));
        end
    end
    ivm_out('  %s %s %s %s %s\n', ivm_pad_left(hdr{1}, wc(1)), ivm_pad_left(hdr{2}, wc(2)), ...
        ivm_pad_left(hdr{3}, wc(3)), ivm_pad_left(hdr{4}, wc(4)), ivm_pad_left(hdr{5}, wc(5)));
    for i = 1:numel(r.kappa_adj_table)
        row = r.kappa_adj_table(i);
        pct = row.kappa_frac * 100.0;
        if using_boot
            ci_lo_use = row.theta_adj - z_crit * se_theta_use;
            ci_hi_use = row.theta_adj + z_crit * se_theta_use;
        else
            ci_lo_use = row.ci_lo;
            ci_hi_use = row.ci_hi;
        end
        ivm_out('  %s %*.4f %*.4f %*.4f %*.4f\n', ivm_pad_left(sprintf('%5.0f%%', pct), wc(1)), ...
            wc(2), row.kappa + 0, wc(3), row.theta_adj, wc(4), ci_lo_use, wc(5), ci_hi_use);
    end
    % A conditional first stage below the OP critical value: c_TM, and so the
    % kappa quantities, are estimated imprecisely (App. F.4 step 3).
    if isequal(ivm_getfield(r, 'F_MZ_T_above_crit', []), false)
        % the labels of the first-stage and kappa lines above (ASCII mode
        % spells the symbols out: tau, theta)
        f_lab = 'F_{M, Z|T, X}';
        names = {TAU, 'c_TM', TH};
        if r.n_instruments > 1
            f_lab = [f_lab, ' (OP eff.)'];
        end
        ivm_out('  Note: %s = %.2f is below the OP critical value (%.2f, %s = %s), so %s is estimated\n', ...
            f_lab, r.F_MZ_T_op, r.F_MZ_T_op_critval, names{1}, ivm_pyrepr(r.op_tau), names{2});
        if ar_shown
            ivm_out('  imprecisely; the AR set for %s above is robust to weak instruments.\n', names{3});
        else
            ivm_out('  imprecisely.\n');
        end
    end
end

% -- R^2-based sensitivity, on request (r2_sensitivity); the "=" signs line
%    up (display width) --
if r2_sensitivity
    ivm_out('\n  R%s-BASED SENSITIVITY (reduced-form robustness)\n', SUP2);
    if r.n_instruments > 1
        ivm_out('  (t-statistics computed under %s SEs; with K > 1, RV/XRV use the first instrument)\n', r.se_type);
    else
        ivm_out('  (t-statistics computed under %s SEs)\n', r.se_type);
    end
    rlab = {sprintf('F (reduced form, = t%s)', SUP2), sprintf('RV(%s)', TH), ...
            sprintf('XRV(%s)', TH)};
    wr = 0;
    for i = 1:3
        wr = max(wr, ivm_dwidth(rlab{i}));
    end
    ivm_out('  %s = %.1f\n', ivm_pad_right_dw(rlab{1}, wr), r.F_reduced_form);
    ivm_out('  %s = %.4f\n', ivm_pad_right_dw(rlab{2}, wr), r.RV_theta);
    ivm_out('  %s = %.4f\n', ivm_pad_right_dw(rlab{3}, wr), r.XRV_theta);
    ivm_out('    (RV_rf=%.4f, RV_fs=%.4f)\n', r.rv_rf, r.rv_fs);
end

if r.lasso && r.lasso_show_selection && ~isempty(r.lasso_info)
    info = r.lasso_info;
    details = info.details;
    ivm_out('\n  LASSO SELECTION DETAILS\n');
    ivm_out('  Penalty method: %s\n', r.lasso_penalty);
    if ~isempty(r.lasso_seed)
        ivm_out('  CV seed: %d\n', r.lasso_seed);
    end
    ivm_print_lasso_set('Y', details.Y, LAM);
    ivm_print_lasso_set('T', details.T, LAM);
    ivm_print_lasso_set('M', details.M, LAM);
    for k = 1:numel(details.Z)
        ivm_print_lasso_set(sprintf('Z[%d]', k - 1), details.Z{k}, LAM);
    end
    ivm_out('  Union S = S_Y %s S_T %s S_M %s S_Z: %d of %d candidates\n', ...
        CUP, CUP, CUP, info.n_selected_union, info.n_total);
    if ~isempty(info.union_names)
        shown = info.union_names(1:min(20, numel(info.union_names)));
        more = info.n_selected_union - numel(shown);
        line = strjoin(shown, ', ');
        if more > 0
            line = [line, sprintf(', ... (+%d more)', more)];
        end
        ivm_out('    %s\n', line);
    end
end

ivm_out('\n%s\n\n', repmat('=', 1, 70));
end

function v = ivm_getfield(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end

function s = ivm_op_str(crit, keff, op_tau, op_alpha)
if isempty(crit) || ~isfinite(crit)
    s = '  |  OP critical value: unavailable';
    return;
end
s = sprintf('  |  OP critical value = %.2f  (K_eff=%.2f, %s=%s, %s=%s)', ...
    crit, keff, ivm_uch(964), ivm_pyrepr(op_tau), ivm_uch(945), ivm_pyrepr(op_alpha));
end

function ivm_print_lasso_set(label, d, LAM)
if isempty(d)
    return;
end
k_it = d.n_iter;
if k_it ~= 1
    pl = 's';
else
    pl = '';
end
if d.converged
    conv = sprintf('converged in %d iteration%s', k_it, pl);
else
    conv = sprintf('did not converge after %d iteration%s', k_it, pl);
end
ivm_out('  S_%s  (%s selected, %s = %.3f, %s)\n', label, ...
    ivm_pad_left(sprintf('%d', d.n_selected), 3), LAM, d.lambda_value, conv);
if d.n_selected > 0
    names = d.names(1:min(20, numel(d.names)));
    more = d.n_selected - numel(names);
    line = strjoin(names, ', ');
    if more > 0
        line = [line, sprintf(', ... (+%d more)', more)];
    end
    ivm_out('    %s\n', line);
end
end

% ============================================================================
%  Bootstrap inference
%
%  'pairs' resamples observations (valid for cross-sectional iid data);
%  'block' resamples whole clusters (valid under cluster dependence). Since
%  2026-10 the estimation sample is what is resampled: rows with missing
%  values in any used column, zero-weight rows and (block) rows without a
%  cluster id are removed first. Each drawn cluster copy is its own cluster,
%  and fixed effects nested in the clusters get one level per drawn copy
%  (as Stata's idcluster()). The bootstrap is not available with lasso=true
%  (naive bootstrap ignores the selection step).
% ============================================================================

function out = ivm_bootstrap(data, Y, T, M, Z, varargin)
%IVM_BOOTSTRAP Pairs/block bootstrap for the IV mediation pipeline
% (Python ivmediate.bootstrap). Testing hooks: 'boot_indices' is a B x n_est
% (pairs) or B x G (block) matrix of 0-BASED draws (rows of the estimation
% sample, or clusters in order of first appearance in the estimation
% sample); when given, the generator is bypassed. 'return_indices' adds the
% draws (out.indices, 0-based) and the estimation sample's 0-based row
% positions in data (out.rows).
defaults = struct();
defaults.X = [];
defaults.dummies = [];
defaults.fe = [];
defaults.se = 'robust';
defaults.clusterid = [];
defaults.alpha = 0.05;
defaults.weights = [];
defaults.wtype = 'aweight';
defaults.op_tau = 0.10;
defaults.op_alpha = 0.05;
defaults.lat = [];
defaults.lon = [];
defaults.spatial_cutoff = [];
defaults.spatial_kernel = [];
defaults.spatial_metric = [];
defaults.method = [];
defaults.B = 999;
defaults.seed = 20260317;
defaults.ci_type = 'basic';
defaults.return_boot_dist = false;
defaults.bootstrap_diagnostics = false;
defaults.progress = true;
defaults.verbose = false;
defaults.quiet = false;
defaults.ascii = [];
defaults.boot_indices = [];
defaults.return_indices = false;
opts = ivm_parse_opts(defaults, varargin, 'ivmediate_bootstrap');
[old_ascii, wstate] = ivm_enter(opts.ascii);
cleanup = onCleanup(@() ivm_leave(old_ascii, wstate)); %#ok<NASGU>
out = ivm_bootstrap_impl(data, Y, T, M, Z, opts);
end

function out = ivm_bootstrap_impl(data, Y, T, M, Z, opts)
return_boot_dist = ivm_check_bool('return_boot_dist', opts.return_boot_dist);
diagnostics = ivm_check_bool('bootstrap_diagnostics', opts.bootstrap_diagnostics);
progress = ivm_check_bool('progress', opts.progress);
verbose = ivm_check_bool('verbose', opts.verbose);
quiet = ivm_check_bool('quiet', opts.quiet);
return_indices = ivm_check_bool('return_indices', opts.return_indices);
se = ivm_check_choice('se', opts.se, {'classical', 'robust', 'cluster', 'spatial'});
ivm_check_open_unit('alpha', opts.alpha);
alpha = double(opts.alpha);
clusterid = ivm_opt_name(opts.clusterid, 'clusterid');
method = opts.method;
if isempty(method) && ~ischar(method)
    if strcmp(se, 'cluster') && ~isempty(clusterid)
        method = 'block';
    else
        method = 'pairs';
    end
end
method = ivm_check_choice('method', method, {'pairs', 'block'});
if strcmp(method, 'block') && isempty(clusterid)
    error(['method=''block'' (cluster bootstrap) requires clusterid. ' ...
           'Provide clusterid or switch to method=''pairs'' for iid data.']);
end
if strcmp(method, 'pairs') && ~isempty(clusterid) && strcmp(se, 'cluster')
    ivm_note(quiet, sprintf(['  Note: method=''pairs'' was selected but clusterid=''%s'' ' ...
        'is provided. Pairs bootstrap is valid only under cross-sectional ' ...
        'independence; if your data are clustered, use method=''block''.'], clusterid));
end
if strcmp(method, 'block') && ~strcmp(se, 'cluster')
    ivm_note(quiet, sprintf(['  Note: block resampling by ''%s'' with se=''%s'' ' ...
        '(the replicates do not use cluster SEs).'], clusterid, se));
end
ci_type = ivm_check_choice('ci_type', opts.ci_type, {'basic', 'percentile', 'studentized'});
B = opts.B;
if ~(isnumeric(B) && isscalar(B) && isreal(B) && isfinite(B) && B == round(B))
    error('B must be an integer; got %s.', ivm_repr(B));
end
B = double(B);
if B < 50
    error('B should be at least 50 (got %d); 999 is recommended.', B);
end
seed = opts.seed;
if ~(isnumeric(seed) && isscalar(seed) && isreal(seed) && isfinite(seed) && ...
     seed == round(seed) && seed >= 0 && seed < 2^53)
    error('seed must be a non-negative integer; got %s.', ivm_repr(seed));
end
seed = double(seed);

% -- Estimation sample: resample only rows the analysis would use --
prep = ivm_prepare_sample(data, Y, T, M, Z, opts.X, opts.dummies, opts.fe, ...
                          se, clusterid, ivm_opt_name(opts.weights, 'weights'), ...
                          ivm_opt_name(opts.lat, 'lat'), ivm_opt_name(opts.lon, 'lon'));
mask = prep.mask;
cols = prep.used;
if strcmp(method, 'block')
    if ~any(strcmp(ivm_data_names(data), clusterid))
        error('column(s) not found in data: %s.', clusterid);
    end
    if ~any(strcmp(cols, clusterid))
        cols{end+1} = clusterid;
    end
end
df_all = ivm_df_from_user(data, cols);
if strcmp(method, 'block')
    cl_all = ivm_df_get(df_all, clusterid);
    if iscell(cl_all)
        cl_ok = ~cellfun(@(s) isempty(strtrim(s)), cl_all);
    else
        cl_ok = ~isnan(cl_all);
    end
    mask = mask & cl_ok;
end
est_rows = find(mask) - 1;                 % 0-based positions in data
df = ivm_df_rows(df_all, mask);
n = ivm_df_nrows(df);
if n == 0
    error('No observations remaining after dropping missing values.');
end

% -- Arguments forwarded to the analysis: X and dummies as the lists
%    resolved on the original data (a glob cannot pick up helper columns);
%    ar_ci = false (AR sets are test inversions, not bootstrap targets) --
if strcmp(se, 'cluster')
    cl_arg = clusterid;
else
    cl_arg = [];
end
aopts = ivm_analysis_defaults();
aopts.X = prep.X;
aopts.dummies = prep.dummies;
aopts.fe = opts.fe;
aopts.se = se;
aopts.clusterid = cl_arg;
aopts.alpha = alpha;
aopts.weights = opts.weights;
aopts.wtype = opts.wtype;
aopts.op_tau = opts.op_tau;
aopts.op_alpha = opts.op_alpha;
aopts.ar_ci = false;
aopts.lat = opts.lat;
aopts.lon = opts.lon;
aopts.spatial_cutoff = opts.spatial_cutoff;
aopts.spatial_kernel = opts.spatial_kernel;
aopts.spatial_metric = opts.spatial_metric;
aopts.quiet = quiet;

% -- Original-data fit (point estimates, analytic SEs and CIs) --
if verbose
    ivm_out('  Fitting original sample...\n');
end
r_original = ivm_analysis_impl(df, Y, T, M, Z, aopts, false);

params = {'total_hat', 'tau_hat', 'theta_hat', 'gamma_hat', ...
          'indirect_hat', 'kappa_star', 'rho_TY_star'};
if diagnostics
    params = [params, {'RV_theta', 'XRV_theta'}];
end
result_key = struct('total_hat', 'total_hat', 'tau_hat', 'tau_hat', ...
    'theta_hat', 'theta_hat', 'gamma_hat', 'gamma_hat', ...
    'indirect_hat', 'indirect', 'kappa_star', 'kappa_star', ...
    'rho_TY_star', 'rho_TY_star', 'RV_theta', 'RV_theta', 'XRV_theta', 'XRV_theta');

point = struct();
for i = 1:numel(params)
    point.(params{i}) = r_original.(result_key.(params{i}));
end
analytic_se = struct();
analytic_se.theta_hat = r_original.se_theta;
analytic_se.tau_hat = r_original.se_tau;
analytic_se.gamma_hat = r_original.se_gamma;
analytic_se.total_hat = r_original.se_total;
if isfinite(r_original.se_indirect)
    analytic_se.indirect_hat = r_original.se_indirect;
else
    analytic_se.indirect_hat = [];
end
analytic_se.kappa_star = [];
analytic_se.rho_TY_star = [];
if diagnostics
    analytic_se.RV_theta = [];
    analytic_se.XRV_theta = [];
end
analytic_ci = struct();
analytic_ci.theta_hat = r_original.ci_theta;
analytic_ci.tau_hat = r_original.ci_tau;
analytic_ci.gamma_hat = r_original.ci_gamma;
analytic_ci.total_hat = r_original.ci_total;
if all(isfinite(r_original.ci_indirect))
    analytic_ci.indirect_hat = r_original.ci_indirect;
else
    analytic_ci.indirect_hat = [];
end
analytic_ci.kappa_star = [];
analytic_ci.rho_TY_star = [];
if diagnostics
    analytic_ci.RV_theta = [];
    analytic_ci.XRV_theta = [];
end

% Replicates need their own SEs only for studentized CIs or RV/XRV.
need_se = strcmp(ci_type, 'studentized') || diagnostics;

% -- Set up resampling (on the estimation sample) --
G = 0;
if strcmp(method, 'block')
    [cl_codes, G] = ivm_factorize(ivm_df_get(df, clusterid));
    cluster_rows = cell(1, G);
    [~, order] = sort(cl_codes);           % stable: rows keep their order
    counts = accumarray(cl_codes, 1, [G, 1]);
    ends = cumsum(counts);
    starts = ends - counts + 1;
    for g = 1:G
        cluster_rows{g} = order(starts(g):ends(g));
    end
    cluster_sizes = counts(:).';
    % FE nested in the clusters: one level per drawn copy
    nested_fe = {};
    for c = 1:numel(prep.fe_cols)
        fc = prep.fe_cols{c};
        [fe_codes, n_lev] = ivm_factorize(ivm_df_get(df, fc));
        pairs = unique([fe_codes, cl_codes], 'rows');
        if size(pairs, 1) == numel(unique(pairs(:, 1)))
            nested_fe(end+1, :) = {fc, fe_codes, n_lev}; %#ok<AGROW>
        end
    end
end

boot_indices = opts.boot_indices;
if ~isempty(boot_indices)
    if strcmp(method, 'pairs')
        expected_cols = n;
    else
        expected_cols = G;
    end
    if size(boot_indices, 1) < B || size(boot_indices, 2) ~= expected_cols
        error('boot_indices must be a %d x %d matrix of 0-based indices.', B, expected_cols);
    end
    bi = double(boot_indices(1:B, :));
    if ~isnumeric(boot_indices) || any(bi(:) ~= round(bi(:))) || any(bi(:) < 0) || ...
            any(bi(:) > expected_cols - 1)
        error('boot_indices must contain integers in 0..%d (0-based).', expected_cols - 1);
    end
end
if return_indices
    if strcmp(method, 'pairs')
        idx_record = zeros(B, n);
    else
        idx_record = zeros(B, G);
    end
end

studentized = strcmp(ci_type, 'studentized');
st_params = {'theta_hat', 'tau_hat', 'gamma_hat', 'total_hat', 'indirect_hat'};
boot_dist = struct();
for i = 1:numel(params)
    boot_dist.(params{i}) = NaN(B, 1);
end
if studentized
    boot_se_replicate = struct();
    for i = 1:numel(st_params)
        boot_se_replicate.(st_params{i}) = NaN(B, 1);
    end
end
n_failed = 0;
if progress
    progress_step = max(floor(B / 20), 1);
else
    progress_step = B + 1;
end
if progress && ~quiet
    ivm_out('  Bootstrap (%s, B=%d)...\n', method, B);
end
aopts_b = aopts;
aopts_b.quiet = true;
if strcmp(method, 'block') && strcmp(se, 'cluster')
    aopts_b.clusterid = '_boot_cluster';
end

for b = 1:B
    if strcmp(method, 'pairs')
        if ~isempty(boot_indices)
            idx_b = double(boot_indices(b, :)).' + 1;   % 0-based -> 1-based
        else
            idx_b = ivm_prng_uniform_int(seed, b, n, n);
        end
        if return_indices
            idx_record(b, :) = idx_b.' - 1;
        end
        data_b = ivm_df_rows(df, idx_b);
    else
        if ~isempty(boot_indices)
            sel = double(boot_indices(b, :)).' + 1;     % 0-based -> 1-based
        else
            sel = ivm_prng_uniform_int(seed, b, G, G);
        end
        if return_indices
            idx_record(b, :) = sel.' - 1;
        end
        concat_idx = vertcat(cluster_rows{sel});
        data_b = ivm_df_rows(df, concat_idx);
        copy_id = repelem((0:G-1).', cluster_sizes(sel).');
        data_b = ivm_df_set(data_b, '_boot_cluster', copy_id);
        for c = 1:size(nested_fe, 1)
            codes_c = nested_fe{c, 2};
            data_b = ivm_df_set(data_b, nested_fe{c, 1}, ...
                copy_id * nested_fe{c, 3} + (codes_c(concat_idx) - 1));
        end
    end

    ok = true;
    try
        rb = ivm_analysis_impl(data_b, Y, T, M, Z, aopts_b, ~need_se);
    catch err
        ok = false;
        n_failed = n_failed + 1;
        if verbose
            ivm_out('    Replicate %d: failed (%s)\n', b - 1, err.message);
        end
    end
    if ok
        for i = 1:numel(params)
            boot_dist.(params{i})(b) = rb.(result_key.(params{i}));
        end
        if studentized
            boot_se_replicate.theta_hat(b) = rb.se_theta;
            boot_se_replicate.tau_hat(b) = rb.se_tau;
            boot_se_replicate.gamma_hat(b) = rb.se_gamma;
            boot_se_replicate.total_hat(b) = rb.se_total;
            boot_se_replicate.indirect_hat(b) = rb.se_indirect;
        end
    end
    if progress && ~quiet && mod(b, progress_step) == 0
        ivm_out('    %5d/%d  (%4.1f%%, %d successful)\n', b, B, ...
            100 * b / B, b - n_failed);
    end
end

if n_failed >= B
    error(['All %d bootstrap replicates failed. Check that the model ' ...
           'specification is consistent with the data.'], B);
end
if n_failed > 0.10 * B
    ivm_note(quiet, sprintf(['  Warning: %d/%d bootstrap replicates failed (%.1f%%). ' ...
        'Bootstrap CIs may be unreliable.'], n_failed, B, 100 * n_failed / B));
end

% -- Aggregate to bootstrap SE and CI --
boot_se = struct();
boot_ci = struct();
for i = 1:numel(params)
    p = params{i};
    bd = boot_dist.(p);
    bootv = bd(~isnan(bd));
    if numel(bootv) > 1
        boot_se.(p) = std(bootv);             % ddof = 1, like np.std(ddof=1)
    else
        boot_se.(p) = NaN;
    end
    if strcmp(ci_type, 'percentile')
        q = ivm_quantile_linear(bootv, [alpha / 2, 1 - alpha / 2]);
        boot_ci.(p) = [q(1), q(2)];
    elseif strcmp(ci_type, 'basic')
        q = ivm_quantile_linear(bootv, [alpha / 2, 1 - alpha / 2]);
        boot_ci.(p) = [2 * point.(p) - q(2), 2 * point.(p) - q(1)];
    else   % studentized
        if any(strcmp(st_params, p))
            seb = boot_se_replicate.(p);
            ok_b = (~isnan(bd)) & (~isnan(seb));
            if sum(ok_b) < 10
                boot_ci.(p) = [NaN, NaN];
                continue;
            end
            t_stats = (bd(ok_b) - point.(p)) ./ seb(ok_b);
            t_stats = t_stats(isfinite(t_stats));
            if numel(t_stats) < 10
                boot_ci.(p) = [NaN, NaN];
                continue;
            end
            q = ivm_quantile_linear(t_stats, [alpha / 2, 1 - alpha / 2]);
            point_se = analytic_se.(p);
            if isempty(point_se) || ~isfinite(point_se)
                boot_ci.(p) = [NaN, NaN];
            else
                boot_ci.(p) = [point.(p) - q(2) * point_se, ...
                               point.(p) - q(1) * point_se];
            end
        else
            q = ivm_quantile_linear(bootv, [alpha / 2, 1 - alpha / 2]);
            boot_ci.(p) = [2 * point.(p) - q(2), 2 * point.(p) - q(1)];
        end
    end
end

out = struct();
out.point = point;
out.boot_se = boot_se;
out.boot_ci = boot_ci;
out.analytic_se = analytic_se;
out.analytic_ci = analytic_ci;
out.n_failed = n_failed;
out.n_est = n;
out.method = method;
out.B = B;
out.seed = seed;
out.ci_type = ci_type;
out.alpha = alpha;
if return_boot_dist
    out.boot_dist = boot_dist;
end
if return_indices
    out.indices = idx_record;
    out.rows = est_rows;
end
end

function ivm_print_bootstrap(boot, varargin)
%IVM_PRINT_BOOTSTRAP Compact bootstrap-results table (Python parity).
defaults = struct();
defaults.title = 'Bootstrap Mediation';
defaults.stars = true;
defaults.star_cutoffs = [0.01, 0.05, 0.10];
defaults.ascii = [];
popts = ivm_parse_opts(defaults, varargin, 'ivmediate_print_bootstrap');
[old_ascii, wstate] = ivm_enter(popts.ascii);
cleanup = onCleanup(@() ivm_leave(old_ascii, wstate)); %#ok<NASGU>
title_str = ivm_norm_text_scalar(popts.title, 'title');
stars = ivm_check_bool('stars', popts.stars);
star_cutoffs = popts.star_cutoffs;
if stars
    star_cutoffs = ivm_validate_star_cutoffs(star_cutoffs);
end
coef_keys = {'theta_hat', 'tau_hat', 'gamma_hat', 'total_hat', 'indirect_hat'};

TH = ivm_uch(952); GA = ivm_uch(947); TAU = ivm_uch(964); KAP = ivm_uch(954);
RHO = ivm_uch(961); ARR = ivm_uch(8594);
HBAR = ivm_uch(9472); ALPH = ivm_uch(945);

ivm_out('\n%s\n', repmat('=', 1, 70));
ivm_out('  %s  (method=%s, B=%d, seed=%d)\n', title_str, boot.method, ...
    boot.B, boot.seed);
ivm_out('%s\n', repmat('=', 1, 70));
ivm_out('  Bootstrap CI type: %s  (%s, %s = %s)\n', boot.ci_type, ...
    ivm_level_label(boot.alpha), ALPH, ivm_pyrepr(boot.alpha));
if boot.n_failed > 0
    failed_str = sprintf(', %d failed', boot.n_failed);
else
    failed_str = '';
end
ivm_out('  %d successful replicates%s.\n', boot.B - boot.n_failed, failed_str);

if ivm_ascii_state()
    TG = 'theta*gamma';
else
    TG = [TH, GA];
end
label_map = struct();
label_map.theta_hat = sprintf('%s (mediator%sY)', TH, ARR);
label_map.tau_hat = sprintf('%s (direct T%sY)', TAU, ARR);
label_map.gamma_hat = sprintf('%s (T%smediator)', GA, ARR);
label_map.total_hat = sprintf('%s_total', TAU);
label_map.indirect_hat = sprintf('%s (indirect)', TG);
label_map.kappa_star = sprintf('%s*', KAP);
label_map.rho_TY_star = sprintf('%s*_TY', RHO);
label_map.RV_theta = sprintf('RV(%s)', TH);
label_map.XRV_theta = sprintf('XRV(%s)', TH);
pnames = fieldnames(boot.point);
labels = cell(1, numel(pnames));
for i = 1:numel(pnames)
    if isfield(label_map, pnames{i})
        labels{i} = label_map.(pnames{i});
    else
        labels{i} = pnames{i};
    end
end
% labels padded to 18 as in Python; ASCII mode: to the longest ASCII label
wl = ivm_label_width(labels, 18);
ivm_out('  %s\n', repmat(HBAR, 1, 68 + wl - 18));
ivm_out('  %s %s %s %s\n', ivm_pad_right('Parameter', wl), ...
    ivm_pad_left('Estimate', 10), ivm_pad_left('Boot SE', 10), ...
    ivm_pad_left('Boot CI', 26));
for i = 1:numel(pnames)
    p = pnames{i};
    label = labels{i};
    est = boot.point.(p);
    sdv = boot.boot_se.(p);
    ci = boot.boot_ci.(p);
    ci_str = sprintf('[%10.4f, %10.4f]', ci(1), ci(2));
    if stars && any(strcmp(coef_keys, p))
        star_str = ivm_sig_stars(ivm_wald_pvalue(est, sdv), star_cutoffs);
    else
        star_str = '';
    end
    if ~isempty(star_str)
        suffix = [' ', star_str];
    else
        suffix = '';
    end
    ivm_out('  %s %10.4f %10.4f %s%s\n', ivm_pad_right(label, wl), ...
        est, sdv, ivm_pad_left(ci_str, 26), suffix);
end
if stars
    ivm_out('  Sig. codes: *** p<%s, ** p<%s, * p<%s\n', ...
        ivm_format_cut(star_cutoffs(1)), ivm_format_cut(star_cutoffs(2)), ...
        ivm_format_cut(star_cutoffs(3)));
end
ivm_out('%s\n\n', repmat('=', 1, 70));
end
