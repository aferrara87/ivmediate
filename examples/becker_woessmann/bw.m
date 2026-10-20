% Becker and Woessmann (2009): IV mediation estimates (paper, Table 6, column 1).
%
% Run from MATLAB:  run('bw.m')   or in batch mode:  matlab -batch "run('bw.m')"
% (in batch mode on Windows the report uses ASCII symbol names; option 'ascii')
% Uses the ivmediate package; when it is not on the path, the repository's
% matlab/ folder is added.

here = fileparts(mfilename('fullpath'));
if exist('ivmediate', 'file') ~= 2
    addpath(fullfile(here, '..', '..', 'matlab'));
end

% The CSV codes missing values as NA; ivmediate_read_csv reads NA as missing.
data = ivmediate_read_csv(fullfile(here, 'becker_woessmann_2009.csv'));
r = ivmediate(data, 'inctaxpc', 'f_prot', 'f_rw', {'kmwittenberg'}, ...
    'X', 'f_young f_jew f_fem f_ortsgeb f_pruss f_blind f_deaf f_dumb hhsize lnpop gpop', ...
    'se', 'robust');
ivmediate_print(r, 'title', 'Becker and Woessmann (2009)');
