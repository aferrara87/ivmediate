% RUN_MATLAB  MATLAB runner for the cross-language validation harness.
%
% MATLAB counterpart of run_py.py and run_r.R. Runs every config of
% configs.json through the MATLAB package (matlab/, version 2.0.0) on
% data_<data>.csv, read with ivmediate_read_csv:
%   configs        ivmediate + ivmediate_print (smoke test with the
%                  optional kappa table and R^2 section on, Addendum A.2)
%   error_configs  degenerate inputs that must stop with an error of class
%                  "expect"
%   boot_configs   ivmediate_bootstrap + ivmediate_print_bootstrap,
%                  replaying the Python draws
%                  reference_2026-10/bootidx_<name>.csv (rows of the
%                  estimation sample, or clusters in order of first appearance
%                  in it); the estimation sample is checked against
%                  bootidx_<name>_rows.csv.
% The version and the internal helpers (OP table, critical values) come from
% the package's testing hook ivmediate('__internal__').
% A config may name Y, T and M (default y, tvar, mvar); every other field
% except name, data and expect is passed as a name/value option.
%
% Each config's output is a flat object of strings with the keys and the
% formatting of run_py.py: numbers '%.17g', NaN, NULL ([] plays Python's
% None), TRUE/FALSE, lists as key.1, key.2, ..., text joined with '|'.
% Metadata keys start with "_": _status ("ok" or "error"), _error_class (the
% first class in configs.json "error_classes" with a phrase contained in the
% lower-cased message, else "other"), _ERROR ("<identifier>: <message>"),
% _stdout_lines (printed lines), and for bootstrap configs _rows_match.
%
% Run from any folder (paths are resolved from this file's location):
%   matlab -batch "run('<repo>/validation/run_matlab.m')"
%   python -B validation/compare.py matlab
% Optional workspace settings before run(...):
%   RUN_MATLAB_OUT   output path (default: validation/out_matlab.json)
%   RUN_MATLAB_ONLY  cellstr of config names to run (default: all)

rm_here = fileparts(mfilename('fullpath'));
if ~exist('RUN_MATLAB_OUT', 'var') || isempty(RUN_MATLAB_OUT)
    RUN_MATLAB_OUT = fullfile(rm_here, 'out_matlab.json');
end
if ~exist('RUN_MATLAB_ONLY', 'var')
    RUN_MATLAB_ONLY = {};
end
rm_main(rm_here, RUN_MATLAB_OUT, RUN_MATLAB_ONLY);


% ============================================================================
function rm_main(val_dir, out_path, only)
root = fileparts(val_dir);
addpath(fullfile(root, 'matlab'));
ivm = ivmediate('__internal__');   % testing hook: the version and the internal helpers
ref_dir = fullfile(val_dir, 'reference_2026-10');

fprintf(1, '[matlab] %s, ivmediate %s\n', version, ivm.version);

cfg_txt = fileread(fullfile(val_dir, 'configs.json'));
CFG = jsondecode(cfg_txt);
classes = rm_error_classes(CFG.error_classes);

datasets = struct();
dnames = {'main', 'weak', 'ext', 'negvar', 'chain'};
for i = 1:numel(dnames)
    t0 = tic;
    datasets.(dnames{i}) = ivmediate_read_csv(fullfile(val_dir, ['data_', dnames{i}, '.csv']));
    fprintf(1, '[matlab] read data_%s.csv (%.2f s)\n', dnames{i}, toc(t0));
end

ALL = containers.Map('KeyType', 'char', 'ValueType', 'any');
n_unexpected = 0;

% ---- static: version, OP table and computed critical values ----
st = containers.Map('KeyType', 'char', 'ValueType', 'any');
st('version') = ivm.version;
[tab, taus] = ivm.op_table();
for K = 1:size(tab, 1)
    for j = 1:numel(taus)
        st(sprintf('op.%d.%.2f', K, taus(j))) = rm_fmt(tab(K, j));
    end
end
cases = [1, 0.10, 0.05; 5, 0.30, 0.05; 30, 0.05, 0.05; ...
         31, 0.10, 0.05; 1, 0.15, 0.05; 1, 0.10, 0.01];
for i = 1:size(cases, 1)
    key = sprintf('opcv.%d.%s.%s', cases(i, 1), rm_pyfloat(cases(i, 2)), ...
                  rm_pyfloat(cases(i, 3)));
    st(key) = rm_fmt(ivm.op_critical_value(cases(i, 1), cases(i, 2), cases(i, 3)));
end
ALL('_static') = st;

% ---- analysis configs and degenerate inputs ----
configs = [rm_as_cell(CFG.configs), rm_as_cell(CFG.error_configs)];
t_all = tic;
for c = 1:numel(configs)
    cfg = configs{c};
    name = cfg.name;
    if ~isempty(only) && ~any(strcmp(only, name))
        continue;
    end
    data = datasets.(cfg.data);
    [roles, Zc, args, expect] = rm_cfg_args(cfg, cfg_txt);
    t0 = tic;
    r = []; %#ok<NASGU>  (assigned inside evalc)
    try
        txt = evalc(['r = ivmediate(data, roles{1}, roles{2}, roles{3}, Zc, args{:}); ' ...
                     'ivmediate_print(r, ''title'', name, ''kappa_table'', true, ''r2_sensitivity'', true);']);
        out = rm_flatten_result(r);
        out('_status') = 'ok';
        out('_stdout_lines') = sprintf('%d', rm_count_lines(txt));
        msg = 'ok';
    catch err
        out = rm_error_entry(err, classes);
        msg = sprintf('ERROR [%s] %s', out('_error_class'), out('_ERROR'));
    end
    got = 'ok';
    if isKey(out, '_error_class')
        got = out('_error_class');
    end
    if ~strcmp(got, expect)
        n_unexpected = n_unexpected + 1;
        msg = sprintf('%s   <-- UNEXPECTED (expected %s)', msg, expect);
    end
    ALL(name) = out;
    fprintf(1, '[matlab] %s: %s (%.2f s)\n', name, msg, toc(t0));
end

% ---- bootstrap configs (Python draws replayed) ----
bconfigs = rm_as_cell(CFG.boot_configs);
for c = 1:numel(bconfigs)
    bcfg = bconfigs{c};
    name = bcfg.name;
    if ~isempty(only) && ~any(strcmp(only, name))
        continue;
    end
    data = datasets.(bcfg.data);
    idx = rm_read_idx(fullfile(ref_dir, sprintf('bootidx_%s.csv', name)));
    rows_py = rm_read_idx(fullfile(ref_dir, sprintf('bootidx_%s_rows.csv', name)));
    [roles, Zc, args, expect] = rm_cfg_args(bcfg, cfg_txt);
    t0 = tic;
    bt = []; %#ok<NASGU>  (assigned inside evalc)
    try
        txt = evalc(['bt = ivmediate_bootstrap(data, roles{1}, roles{2}, roles{3}, Zc, args{:}, ' ...
                     '''return_boot_dist'', true, ''progress'', false, ''boot_indices'', idx, ' ...
                     '''return_indices'', true); ' ...
                     'ivmediate_print_bootstrap(bt, ''title'', name);']);
        out = rm_flatten_boot(bt);
        out('_status') = 'ok';
        out('_stdout_lines') = sprintf('%d', rm_count_lines(txt));
        rows_ok = isequal(double(bt.rows(:)), double(rows_py(:))) && ...
                  isequal(double(bt.indices), double(idx));
        out('_rows_match') = rm_fmt(rows_ok);
        msg = sprintf('ok (n_est=%d, n_failed=%d, rows match: %d)', bt.n_est, bt.n_failed, rows_ok);
        if ~rows_ok
            n_unexpected = n_unexpected + 1;
            msg = [msg, '   <-- ESTIMATION SAMPLE DIFFERS FROM PYTHON']; %#ok<AGROW>
        end
    catch err
        out = rm_error_entry(err, classes);
        msg = sprintf('ERROR [%s] %s', out('_error_class'), out('_ERROR'));
    end
    got = 'ok';
    if isKey(out, '_error_class')
        got = out('_error_class');
    end
    if ~strcmp(got, expect)
        n_unexpected = n_unexpected + 1;
        msg = sprintf('%s   <-- UNEXPECTED (expected %s)', msg, expect);
    end
    ALL(name) = out;
    fprintf(1, '[matlab] %s: %s (%.2f s)\n', name, msg, toc(t0));
end
fprintf(1, '[matlab] all configs: %.1f s\n', toc(t_all));

rm_write_json(out_path, ALL);
fprintf(1, '[matlab] wrote %s (%d configs; %d with an unexpected status)\n', ...
        out_path, ALL.Count - 1, n_unexpected);
end


% ============================================================================
%  Configs, error classes, indices
% ============================================================================

function c = rm_as_cell(x)
if isempty(x)
    c = {};
elseif iscell(x)
    c = x(:).';
else
    c = num2cell(x(:).');
end
end

function classes = rm_error_classes(ec)
% configs.json "error_classes": [[class, [phrases...]], ...] -> N x 2 cell.
ec = rm_as_cell(ec);
classes = cell(numel(ec), 2);
for i = 1:numel(ec)
    e = ec{i};
    classes{i, 1} = e{1};
    ph = e{2};
    if ischar(ph)
        ph = {ph};
    end
    classes{i, 2} = ph(:).';
end
end

function cls = rm_error_class(msg, classes)
low = lower(msg);
for i = 1:size(classes, 1)
    ph = classes{i, 2};
    for j = 1:numel(ph)
        if contains(low, ph{j})
            cls = classes{i, 1};
            return;
        end
    end
end
cls = 'other';
end

function out = rm_error_entry(err, classes)
out = containers.Map('KeyType', 'char', 'ValueType', 'any');
out('_status') = 'error';
out('_error_class') = rm_error_class(err.message, classes);
ident = err.identifier;
if isempty(ident)
    ident = 'MATLAB:error';
end
out('_ERROR') = sprintf('%s: %s', ident, err.message);
end

function [roles, Zc, args, expect] = rm_cfg_args(cfg, cfg_txt)
% Roles (Y, T, M; defaults y, tvar, mvar), the positional Z and the
% name/value pairs from one config object (name, data, expect removed).
roles = {'y', 'tvar', 'mvar'};
rn = {'Y', 'T', 'M'};
for i = 1:3
    if isfield(cfg, rn{i})
        roles{i} = cfg.(rn{i});
    end
end
Zc = cfg.Z;
if ischar(Zc)
    Zc = {Zc};
end
Zc = Zc(:).';
expect = 'ok';
if isfield(cfg, 'expect')
    expect = cfg.expect;
end
args = {};
fn = fieldnames(cfg);
for i = 1:numel(fn)
    f = fn{i};
    if any(strcmp(f, {'name', 'data', 'Z', 'Y', 'T', 'M', 'expect'}))
        continue;
    end
    v = cfg.(f);
    if iscell(v)
        v = v(:).';
    end
    if strcmp(f, 'expand_override')
        % jsondecode mangles keys such as '_D_dnum_2' into valid field names
        % ('x_D_dnum_2'), so rebuild the override from the raw JSON text as
        % the documented N x 2 cell {name, {ops...}}.
        v = rm_raw_override(cfg_txt, cfg.name);
    end
    args = [args, {f, v}]; %#ok<AGROW>
end
end

function ov = rm_raw_override(cfg_txt, cfg_name)
k = strfind(cfg_txt, sprintf('"name": "%s"', cfg_name));
if isempty(k)
    error('run_matlab: config %s not found in configs.json text.', cfg_name);
end
rest = cfg_txt(k(1):end);
k2 = strfind(rest, '"expand_override"');
rest = rest(k2(1):end);
lb = find(rest == '{', 1);
rb = find(rest == '}', 1);
body = rest(lb + 1:rb - 1);
toks = regexp(body, '"([^"]+)"\s*:\s*\[([^\]]*)\]', 'tokens');
ov = cell(numel(toks), 2);
for i = 1:numel(toks)
    ov{i, 1} = toks{i}{1};
    ops = regexp(toks{i}{2}, '"([^"]+)"', 'tokens');
    ov{i, 2} = cellfun(@(t) t{1}, ops, 'UniformOutput', false);
end
end

function idx = rm_read_idx(path)
% Read a comma-separated integer matrix (one row per line).
txt = strtrim(fileread(path));
lines = regexp(txt, '\r?\n', 'split');
nrow = numel(lines);
vals = sscanf(strrep(txt, ',', ' '), '%d');
if mod(numel(vals), nrow) ~= 0
    error('run_matlab: ragged index file %s.', path);
end
idx = reshape(vals, [], nrow).';
end

function n = rm_count_lines(txt)
% Python str.splitlines() count: a trailing newline does not add a line.
if isempty(txt)
    n = 0;
    return;
end
parts = regexp(txt, '\r\n|\n|\r', 'split');
if isempty(parts{end})
    parts = parts(1:end-1);
end
n = numel(parts);
end


% ============================================================================
%  Flattening (mirrors run_py.py flatten_result / flatten_boot)
% ============================================================================

function s = rm_fmt(v)
% One value -> string, following run_py.py put()/fnum().
if ischar(v)
    s = v;
elseif isa(v, 'string')
    s = char(v);
elseif isempty(v)
    s = 'NULL';                              % [] plays Python's None
elseif numel(v) ~= 1
    s = sprintf('<non-scalar %s %s>', class(v), mat2str(size(v)));
elseif islogical(v)
    if v
        s = 'TRUE';
    else
        s = 'FALSE';
    end
elseif isnumeric(v)
    v = double(v);
    if isnan(v)
        s = 'NaN';
    elseif isinf(v)
        if v > 0
            s = 'inf';                       % Python '%.17g' spelling
        else
            s = '-inf';
        end
    else
        s = sprintf('%.17g', v);
    end
else
    s = sprintf('<unsupported %s>', class(v));
end
end

function s = rm_pyfloat(v)
% Python str() of a short float such as 0.1 or 0.05.
s = sprintf('%g', v);
end

function put(out, key, v)
out(key) = rm_fmt(v); %#ok<NASGU>  (containers.Map is a handle)
end

function put_vec(out, key, v)
% Scalars as `key`, vectors (K >= 2) as `key.1`, `key.2`, ...
if numel(v) <= 1
    put(out, key, v);
else
    for i = 1:numel(v)
        put(out, sprintf('%s.%d', key, i), v(i));
    end
end
end

function v = rm_get(r, k)
% Python r.get(k): a missing field reads as None.
if isfield(r, k)
    v = r.(k);
else
    fprintf(1, '[matlab]   note: result has no field %s\n', k);
    v = [];
end
end

function out = rm_flatten_result(r)
out = containers.Map('KeyType', 'char', 'ValueType', 'any');
keys = {'n', 'n_instruments', 'n_covariates', 'n_dummies', 'se_type', 'alpha', ...
        'n_fe_groups', 'n_dropped_singletons', 'df_fe_absorbed', 'weighted', 'wtype', ...
        'spatial_cutoff', 'spatial_kernel', 'spatial_metric', ...
        'theta_hat', 'se_theta', 'tau_hat', 'se_tau', 'gamma_hat', 'se_gamma', ...
        'total_hat', 'se_total', 'indirect', 'se_indirect', 'cov_gamma_theta', ...
        'decomp_check', 'decomp_error', ...
        'F_TZ', 'F_MZ_T', 'F_TZ_op', 'F_MZ_T_op', ...
        'F_TZ_op_critval', 'F_MZ_T_op_critval', 'op_tau', 'op_alpha', ...
        'c_TM_hat', 'rho_TM_hat', 'sigma_eT', 'sigma_eY', ...
        'kappa_star', 'kappa_star_CI', 'rho_TY_star', ...
        'f_reduced_form', 'f_first_stage_cond', ...
        'RV_theta', 'XRV_theta', 'rv_rf', 'rv_fs', 'xrv_rf', 'xrv_fs', ...
        'lasso', 'lasso_penalty', ...
        'n_clusters', 'K_eff_TZ', 'K_eff_MZ_T', ...
        'kappa_star_CI_informative', 'F_MZ_T_above_crit'};
for i = 1:numel(keys)
    put(out, keys{i}, rm_get(r, keys{i}));
end
put_vec(out, 'pi1', rm_get(r, 'pi1'));
put_vec(out, 'se_pi1', rm_get(r, 'se_pi1'));
nested = rm_get(r, 'fe_nested_in_cluster');
if isempty(nested)
    put(out, 'fe_nested_in_cluster', 'NULL');
else
    tf = {'FALSE', 'TRUE'};
    put(out, 'fe_nested_in_cluster', strjoin(tf(double(logical(nested)) + 1), ','));
end
if isempty(r.fe)
    put(out, 'fe', 'NULL');
else
    fe = r.fe;
    if ischar(fe)
        fe = {fe};
    end
    put(out, 'fe', strjoin(fe, ','));
end
cinames = {'ci_theta', 'ci_tau', 'ci_gamma', 'ci_total', 'ci_indirect'};
for i = 1:numel(cinames)
    ci = r.(cinames{i});
    put(out, [cinames{i}, '.1'], ci(1));
    put(out, [cinames{i}, '.2'], ci(2));
end
for i = 1:numel(r.F_MZ_T_individual)
    put(out, sprintf('F_MZ_T_individual.%d', i), r.F_MZ_T_individual(i));
end
kk = {'kappa', 'kappa_frac', 'theta_adj', 'ci_lo', 'ci_hi'};
for i = 1:numel(r.kappa_adj_table)
    row = r.kappa_adj_table(i);
    for j = 1:numel(kk)
        put(out, sprintf('kappa_adj.%d.%s', i, kk{j}), row.(kk{j}));
    end
end
arn = {'ar_theta', 'ar_gamma', 'ar_total'};
for i = 1:numel(arn)
    rm_flat_ar(out, arn{i}, r.(arn{i}));
end
epn = {'eps_T_hat', 'eps_M_hat', 'eps_Y_hat'};
for i = 1:numel(epn)
    rm_flat_eps(out, epn{i}, r.(epn{i}));
end
li = r.lasso_info;
if ~isempty(li)
    put(out, 'lasso_info.n_selected_union', li.n_selected_union);
    put(out, 'lasso_info.n_total', li.n_total);
    put(out, 'lasso_info.union_names', strjoin(li.union_names, '|'));
    det = li.details;
    rm_flat_lasso_detail(out, 'lasso_info.Y', det.Y);
    rm_flat_lasso_detail(out, 'lasso_info.T', det.T);
    rm_flat_lasso_detail(out, 'lasso_info.M', det.M);
    for k = 1:numel(det.Z)
        rm_flat_lasso_detail(out, sprintf('lasso_info.Z%d', k), det.Z{k});
    end
end
end

function rm_flat_ar(out, pre, ar)
if isempty(ar)
    put(out, pre, 'NULL');
    return;
end
kk = {'shape', 'method', 'contains_zero', 'cv', 'df_chi2', 'alpha', ...
      'stat_at_zero', 'pval_at_zero', 'lower', 'upper', ...
      'left_upper', 'right_lower'};
for j = 1:numel(kk)
    put(out, [pre, '.', kk{j}], ar.(kk{j}));
end
% Python stores segments only for the grid path (None for Fieller shapes).
if isfield(ar, 'has_segments') && ar.has_segments
    segs = ar.segments;
    put(out, [pre, '.n_segments'], numel(segs));
    for i = 1:numel(segs)
        put(out, sprintf('%s.seg%d.lo', pre, i), segs(i).lo);
        put(out, sprintf('%s.seg%d.hi', pre, i), segs(i).hi);
        put(out, sprintf('%s.seg%d.open_left', pre, i), logical(segs(i).open_left));
        put(out, sprintf('%s.seg%d.open_right', pre, i), logical(segs(i).open_right));
    end
end
end

function rm_flat_eps(out, pre, e)
e = double(e(:));
put(out, [pre, '.mean'], mean(e));
put(out, [pre, '.sd'], sqrt(mean((e - mean(e)).^2)));
for i = 1:3
    put(out, sprintf('%s.v%d', pre, i), e(i));
end
end

function rm_flat_lasso_detail(out, pre, d)
if isempty(d)
    put(out, pre, 'NULL');
    return;
end
put(out, [pre, '.n_selected'], d.n_selected);
put(out, [pre, '.lambda_value'], d.lambda_value);
put(out, [pre, '.converged'], logical(d.converged));
put(out, [pre, '.n_iter'], d.n_iter);
put(out, [pre, '.method'], d.method);
put(out, [pre, '.names'], strjoin(d.names, '|'));
cs = cell(1, numel(d.coefs));
for i = 1:numel(d.coefs)
    cs{i} = sprintf('%.10g', d.coefs(i));
end
put(out, [pre, '.coefs'], strjoin(cs, '|'));
is = cell(1, numel(d.indices));
for i = 1:numel(d.indices)
    is{i} = sprintf('%d', d.indices(i));    % already 1-based in MATLAB
end
put(out, [pre, '.indices'], strjoin(is, '|'));
end

function out = rm_flatten_boot(bt)
out = containers.Map('KeyType', 'char', 'ValueType', 'any');
kk = {'n_failed', 'method', 'B', 'ci_type', 'alpha', 'n_est'};
for i = 1:numel(kk)
    put(out, kk{i}, bt.(kk{i}));
end
pn = fieldnames(bt.point);
for i = 1:numel(pn)
    put(out, ['point.', pn{i}], bt.point.(pn{i}));
end
pn = fieldnames(bt.boot_se);
for i = 1:numel(pn)
    put(out, ['boot_se.', pn{i}], bt.boot_se.(pn{i}));
end
pn = fieldnames(bt.boot_ci);
for i = 1:numel(pn)
    ci = bt.boot_ci.(pn{i});
    put(out, ['boot_ci.', pn{i}, '.1'], ci(1));
    put(out, ['boot_ci.', pn{i}, '.2'], ci(2));
end
pn = fieldnames(bt.analytic_se);
for i = 1:numel(pn)
    put(out, ['analytic_se.', pn{i}], bt.analytic_se.(pn{i}));
end
if isfield(bt, 'boot_dist')
    pn = fieldnames(bt.boot_dist);
    for i = 1:numel(pn)
        vec = bt.boot_dist.(pn{i});
        ok = vec(~isnan(vec));
        put(out, sprintf('boot_dist.%s.mean', pn{i}), mean(ok));
        put(out, sprintf('boot_dist.%s.sd', pn{i}), std(ok));   % ddof = 1
        for j = 1:3
            put(out, sprintf('boot_dist.%s.v%d', pn{i}, j), vec(j));
        end
    end
end
end


% ============================================================================
%  JSON writer (all values are strings; non-ASCII escaped as \uXXXX like
%  Python's json.dump default, so the file is pure ASCII)
% ============================================================================

function rm_write_json(path, ALL)
fid = fopen(path, 'w');
if fid < 0
    error('run_matlab: cannot write %s', path);
end
cleaner = onCleanup(@() fclose(fid));
names = sort(ALL.keys());
fprintf(fid, '{\n');
for i = 1:numel(names)
    m = ALL(names{i});
    ks = sort(m.keys());
    fprintf(fid, '"%s": {\n', rm_jesc(names{i}));
    for j = 1:numel(ks)
        if j < numel(ks)
            sep = ',';
        else
            sep = '';
        end
        fprintf(fid, '"%s": "%s"%s\n', rm_jesc(ks{j}), rm_jesc(m(ks{j})), sep);
    end
    if i < numel(names)
        fprintf(fid, '},\n');
    else
        fprintf(fid, '}\n');
    end
end
fprintf(fid, '}\n');
end

function o = rm_jesc(s)
s = char(s);
d = double(s);
if all(d >= 32 & d <= 126 & s ~= '"' & s ~= '\')
    o = s;
    return;
end
parts = cell(1, numel(s));
for i = 1:numel(s)
    if s(i) == '"'
        parts{i} = '\"';
    elseif s(i) == '\'
        parts{i} = '\\';
    elseif d(i) < 32 || d(i) > 126
        parts{i} = sprintf('\\u%04x', d(i));
    else
        parts{i} = s(i);
    end
end
o = [parts{:}];
end
