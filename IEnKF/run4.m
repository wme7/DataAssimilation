clear;
defaults;
NEXPERIMENT = 4;

% Specify the model here. Set the rest in <model>_setup.m.

model = 2;

% Do not edit below.

if model == 1
    L3_setup;
elseif model == 2
    L40_setup;
end

ncycle = prm.ncycle;
nexperiment = prm.nexperiment;
verbose = prm.verbose;
seed = prm.seed;
nitermax = prm.nitermax;

if ~isfield(prm, 'obs_nstep')
    prm.obs_nstep = prm.cycle_nstep;
end
if ~isfield(prm, 'doplot')
    doplot = 0;
else
    doplot = prm.doplot;
end

prm_print(prm);
fprintf('\n');

if isnan(prm.inflation)
    fprintf('\n  The inflation magnitude is set to NaN in this experiment\n');
    fprintf('  It needs to be set to a meaningful value to proceed\n');
    return
end

if doplot
    f = figure;
    firstplot = 1;
end

stats_all = [];
for experiment = 1 : NEXPERIMENT
    stats = [];
    prm.seed = experiment;
    randn('state', prm.seed);
    rand('state', prm.seed);
    [x_true, E] = prm.generate(prm);
    fprintf('\n  seed = %d:', prm.seed);
    for cycle = 1 : ncycle
        [status, x_true, E, stats_now] = ienkf_cycle(prm, x_true, E);
        if cycle > prm.nspinup
            stats = stats_merge(stats, stats_now);
            if isfield(prm, 'stats_custom')
                stats = prm.stats_custom(stats, x_true, E);
            end
        end
        if verbose
            if status == 0
                fprintf('.');
            elseif status == 1
                fprintf(',');
            else
                fprintf('|');
            end
        elseif mod(cycle, 10000) == 0
            fprintf('|');
        elseif mod(cycle, 1000) == 0
            fprintf('.');
        end
        if status > 1
            fprintf('\n  error: run(): fatal error %d at cycle = %d', status, cycle);
            break; % for cycle
        end
    end
    if ~isempty(stats)
        fprintf('\n');
        stats_print(prm, stats);
        stats_all = stats_merge(stats_all, stats);
    end
end
fprintf('\n  overall:\n');
stats_print(prm, stats_all);
fprintf('\n');
