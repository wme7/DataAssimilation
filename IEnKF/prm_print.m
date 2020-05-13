function [] = prm_print(prm)
    
    fprintf('  model = "%s", method = "%s"', prm.modelname, prm.method);
    if strcmp(prm.method, 'IEKF')
        fprintf(', C = %.3g', prm.C);
    end
    if prm.cycle_nstep ~= 1
        fprintf(', T = %d', prm.cycle_nstep);
    end
    if prm.ncycle > 1
        fprintf(', obs_var = %.3g, nstep = %d - %d', prm.obs_variance, prm.ncycle, prm.nspinup);
    elseif prm.ncycle == 1
        fprintf(', obs_var = %.3g, nexperiment = %d', prm.obs_variance, prm.nexperiment);
    end
    fprintf(', m = %d', prm.m);
    fprintf(', infl = %.3g', prm.inflation);
    % uncomment when calculating for Figure 2
    %fprintf(', infl = %.3g', (prm.inflation - 1) * 12 / prm.cycle_nstep);
    if prm.obs_nstep ~= prm.cycle_nstep
        fprintf(', obs_nstep = %d', prm.obs_nstep);
    end
    if prm.p ~= prm.n
        fprintf(', p = %d', prm.p);
    end
    if isfield(prm, 'rotate') & prm.rotate
        fprintf(', rotate = %d', prm.rotate);
    end
    if ~strcmp(prm.method, 'EnKF')
        fprintf(', nitermax = %d', prm.nitermax);
    end
    if ~strcmp(prm.method, 'EnKF') & isfield(prm, 'epsilon') & prm.epsilon ~= 0
        fprintf(', epsilon = %.3g', prm.epsilon);
    end
    if prm.verbose
        fprintf('\n');
    else
        fprintf(':');
    end
    
    return
    