defaults;

METHOD = 1;
REVERSE = 0;

prm.generate = @L40_generate;
prm.propagate = @L40_step_c;
prm.collapseratio = 0.2;
prm.n = 40;
prm.m = 40;
prm.observe = @L40_observe_1;
prm.p = length(prm.observe([1 : prm.n], 0));
prm.obs_variance = 1;
prm.ncycle = 20500;
prm.nspinup = 500;
prm.nitermax = 30;

switch METHOD
  case 1
    prm.method = 'EnKF';
  case 2
    prm.method = 'IEnKF';
  case 3
    prm.method = 'IEKF';
end

if strcmp(prm.method, 'IEnKF')
    infls = [0.08 0.10 0.10 0.10 0.12 0.12 0.12 0.12 0.14 0.14 0.16 0.18 0.20 0.20 0.24 0.26];
elseif strcmp(prm.method, 'EnKF')
    infls = [0.08 0.14 0.20 0.26 0.32 0.34 0.34 0.34 0.40 0.48 0.48 0.50 0.50 0.52 0.52 0.52];
elseif strcmp(prm.method, 'IEKF')
    infls = [0.34 0.42 0.40 0.40 0.40 0.40 0.40 0.40 0.40 0.42 0.44 0.52 0.62 0.72 0.72 0.72];
end

N = length(infls);

for i = 1 : N
    if ~REVERSE
        infl = infls(i);
        prm.cycle_nstep = i;
    else
        infl = infls(N + 1 - i);
        prm.cycle_nstep = N + 1 - i;
    end;
    prm.obs_nstep = prm.cycle_nstep;

    if ~exist('out', 'dir')
        system('mkdir out');
    end
    fnameout = sprintf('out/%s-%d-%.2f.mat', prm.method, prm.cycle_nstep, infl);
    fprintf('  %s: ', fnameout);
    if exist(fnameout, 'file')
        fprintf('done\n');
        continue;
    end
        
    prm.inflation = 1 + infl * prm.cycle_nstep / 12;

    randn('state', prm.seed);
    rand('state', prm.seed);

    stats = [];
    [x_true, E] = prm.generate(prm);
    for cycle = 1 : prm.ncycle
        [status, x_true, E, stats_now] = ienkf_cycle(prm, x_true, E, cycle);

        if cycle > prm.nspinup
            stats = stats_merge(stats, stats_now);
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
        
        if mod(cycle, 6000) == 0
            fprintf('|%.3g|', mean(stats.rmse_a));
        elseif mod(cycle, 1000) == 0
            fprintf('.');
        end
    end
        
    save(fnameout, 'stats', 'prm');
    fprintf('|%.3g done\n', mean(stats.rmse_a));
end    
