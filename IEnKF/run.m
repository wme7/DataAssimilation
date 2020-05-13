clear;
defaults;

% Specify the model here. Set the rest in <model>_setup.m.

model = 1;

% do not edit below

if model == 1
    L3_setup;
elseif model == 2
    L40_setup;
end

ncycle = prm.ncycle;
verbose = prm.verbose;
seed = prm.seed;
nitermax = prm.nitermax;
if ~isfield(prm, 'obs_nstep')
    prm.obs_nstep = prm.cycle_nstep;
end
if ~isfield(prm, 'doplot')
    doplot = 1;
else
    doplot = prm.doplot;
end

prm_print(prm);

randn('state', prm.seed);
rand('state', prm.seed);

if doplot
    f = figure;
    firstplot = 1;
end

stats = [];
[x_true, E] = prm.generate(prm);
for cycle = 1 : ncycle
    [status, x_true, E, stats_now] = ienkf_cycle(prm, x_true, E, cycle);

    if doplot & stats_now.niter > 0
        if cycle > prm.nspinup
            figure(f);
            if ~firstplot
                hold on;
            end
            col = 'b';
            x1 = stats_now.rmse(1);
            x2 = stats_now.rmse(stats_now.niter);
            if x1 > x2
                plotfunc([cycle cycle], [x1 x2], 'color', [0.75 0.75 0.75]);
            else
                plotfunc([cycle cycle], [x1 x2], 'color', [1 0.75 0.75]);
            end
            if firstplot
                hold on;
            end
            for iter = 1 : stats_now.niter - 1
                plotfunc(cycle, stats_now.rmse(iter), '+', 'color', col, 'markersize', 5);
                col = nextcolor(col);
            end
            plotfunc(cycle, x2, 'k.', 'markersize', 9);
            hold off;
            if firstplot
                xlabel('step');
                ylabel('RMSE');
                title(sprintf('%s.%d', prm.modelname, SETUP));
            end
            if firstplot
                firstplot = 0;
            end
            pause(0.01);
        end
    end
            
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
            fprintf('|%.3g|', mean(stats.rmse_a));
    elseif mod(cycle, 1000) == 0
        fprintf('.');
    end
    if status > 1
        if ~isempty(stats )
            fprintf('\n');
            stats_print(prm, stats);
        end
        diagnosis = '';
        if status == 2
            diagnosis = ' (numerics)';
        elseif status == 3
            diagnosis = ' (divergence)';
        end
        error(sprintf('\n  error: run(): fatal error %d%s at cycle = %d', status, diagnosis, cycle));
    end
end

if prm.calc_norms
    figure;
    plot(stats.normA2_f ./ stats.normA1_f, 'b.-');
    hold on;
    plot(stats.normT, 'r.-');
    ylim([0 40]);
    set(gca,'DataAspectRatio', [12 5 1]);
    legend('|| A_2^0 || / || A_1^0 ||', '|| T^{-1} ||', 'location', 'northwest');
    xlabel('cycle');
    ylabel('norm');
end
fprintf('\n');

if ~isempty(stats)
    stats_print(prm, stats);
else
    fprintf('\n  warning: could not print run stats: empty variable\n');
end
