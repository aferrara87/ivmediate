function tests = test_bw
%TEST_BW  Smoke test: Becker and Woessmann (2009), the paper's Table 6, column 1.
%
%   Plain run, from any folder (an error if a value differs):
%       run('<repo>/matlab/tests/test_bw.m')
%       matlab -batch "run('<repo>/matlab/tests/test_bw.m')"
%   With the MATLAB unit-test framework:
%       runtests('<repo>/matlab/tests')
%
%   The data are the Becker-Woessmann extract in examples/becker_woessmann/
%   of the repository. When ivmediate is not on the path, the package folder
%   (the parent of this folder) is added.

if nargout > 0
    tests = functiontests(localfunctions);
    return;
end
[bad, r] = bw_mismatches();
if ~strcmp(r.version, '2.0.0')
    bad{end+1} = sprintf('version %s (expected 2.0.0)', r.version);
end
if ~isempty(bad)
    error('test_bw:mismatch', 'Becker-Woessmann Table 6, column 1 not reproduced: %s', ...
          strjoin(bad, '; '));
end
fprintf('ivmediate %s: Becker-Woessmann Table 6, column 1 reproduced (%d values)\n', ...
        r.version, size(bw_expected(), 1) + 2);
end

function test_table6_col1(tc)
bad = bw_mismatches();
tc.verifyEmpty(bad, strjoin(bad, '; '));
end

function test_version(tc)
[~, r] = bw_mismatches();
tc.verifyEqual(r.version, '2.0.0');
end

function [bad, r] = bw_mismatches()
% The paper's values that the estimates do not reproduce (empty if none).
here = fileparts(mfilename('fullpath'));
if exist('ivmediate', 'file') ~= 2
    addpath(fullfile(here, '..'));
end
data = ivmediate_read_csv(fullfile(here, '..', '..', 'examples', 'becker_woessmann', ...
                                   'becker_woessmann_2009.csv'));
r = ivmediate(data, 'inctaxpc', 'f_prot', 'f_rw', {'kmwittenberg'}, ...
    'X', 'f_young f_jew f_fem f_ortsgeb f_pruss f_blind f_deaf f_dumb hhsize lnpop gpop', ...
    'quiet', true);
E = bw_expected();
bad = {};
for i = 1:size(E, 1)
    v = double(r.(E{i, 1}));
    if ~(abs(round(v, E{i, 3}) - E{i, 2}) < 1e-9)
        bad{end+1} = sprintf('%s = %.6g (paper %g)', E{i, 1}, v, E{i, 2}); %#ok<AGROW>
    end
end
ar = r.ar_theta;
if ~strcmp(ar.shape, 'interval') || abs(round(ar.lower, 2) - 0.43) > 1e-9 || ...
        abs(round(ar.upper, 2) - 10.20) > 1e-9
    bad{end+1} = sprintf('AR set for theta: %s [%.4f, %.4f] (paper interval [0.43, 10.20])', ...
                         ar.shape, ar.lower, ar.upper);
end
end

function E = bw_expected()
% Table 6, column 1, as printed in the paper (robust SEs, one instrument):
% result field, value, decimals.
E = {'n',            426,  0
     'total_hat',   0.586, 3
     'se_total',    0.242, 3
     'tau_hat',    -0.299, 3
     'se_tau',      0.214, 3
     'theta_hat',   4.709, 3
     'se_theta',    2.164, 3
     'gamma_hat',   0.188, 3
     'se_gamma',    0.027, 3
     'indirect',    0.886, 3
     'se_indirect', 0.408, 3
     'F_TZ',        75.69, 2
     'F_MZ_T',      16.36, 2
     'rho_TY_star', -0.257, 3};
end
