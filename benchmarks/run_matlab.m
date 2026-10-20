function run_matlab(sizes)
%RUN_MATLAB  Speed comparison, MATLAB (benchmarks/README.md).
%
%   Times ivmediate on data/bench_<N>.csv with the defaults and with the
%   post-double-selection LASSO: one untimed warm-up, then 3 timed runs (only
%   the estimation call; loading the data and printing are excluded). Writes
%   results/matlab.csv and results/env_matlab.txt.
%
%   matlab -batch "run_matlab"                (all sizes)
%   matlab -batch "run_matlab(1000)"          (selected sizes)

here = fileparts(mfilename('fullpath'));
if exist('ivmediate', 'file') ~= 2          % a checkout without the package on the path
    addpath(fullfile(here, '..', 'matlab'));
end
user_sizes = nargin > 0;
if ~user_sizes
    sizes = [1000 10000 100000];
end
X = strjoin(arrayfun(@(j) sprintf('x%d', j), 1:50, 'UniformOutput', false), ' ');
specs = {'default', false; 'lasso', true};
reps = 3;

lines = {'language,n,spec,t1,t2,t3,median,theta,gamma,tau,total,se_type'};
for n = sizes
    data = ivmediate_read_csv(fullfile(here, 'data', sprintf('bench_%d.csv', n)));
    for s = 1:size(specs, 1)
        fit(data, X, specs{s, 2});                       % warm-up
        times = zeros(1, reps);
        for i = 1:reps
            t0 = tic;
            r = fit(data, X, specs{s, 2});
            times(i) = toc(t0);
        end
        lines{end+1} = sprintf('MATLAB,%d,%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%s', ...
            n, specs{s, 1}, times, median(times), r.theta_hat, r.gamma_hat, r.tau_hat, ...
            r.total_hat, r.se_type); %#ok<AGROW>
        fprintf('MATLAB n=%6d %-7s median %8.3f s\n', n, specs{s, 1}, median(times));
    end
end

if ~exist(fullfile(here, 'results'), 'dir')
    mkdir(fullfile(here, 'results'));
end
suffix = '';
if user_sizes
    suffix = ['_', strjoin(arrayfun(@num2str, sizes, 'UniformOutput', false), '_')];
end
write_lines(fullfile(here, 'results', ['matlab', suffix, '.csv']), lines);
write_lines(fullfile(here, 'results', 'env_matlab.txt'), ...
    {sprintf('MATLAB %s; ivmediate %s; %s; computational threads %d', version, r.version, ...
             strtrim(version('-blas')), maxNumCompThreads)});
end

function r = fit(data, X, lasso)
evalc(['r = ivmediate(data, ''y'', ''t'', ''m'', {''z''}, ''X'', X, ''fe'', ''g'', ', ...
       '''lasso'', lasso);']);
end

function write_lines(path, lines)
fid = fopen(path, 'w');
fprintf(fid, '%s\n', lines{:});
fclose(fid);
end
