prm.modelname = 'L3';
prm.seed = 2;
if 0 == 1
    prm.propagate = @L3_step;
else
    prm.propagate = @L3_step_c; % a faster option
end
prm.observe = @L3_observe;
prm.generate = @L3_generate;
prm.n = 3;
prm.p = 3;
plotfunc = @semilogy;
prm.collapseratio = 0.005;

% set the experiment # here

SETUP = 2;

% set the output to your convenience

prm.doplot = 0;
prm.verbose = 0;
prm.calc_norms = 0;

% do not edit below

if SETUP >= 1 & SETUP <= 8
    prm.cycle_nstep = 25;
    prm.obs_variance = 2;
    prm.nspinup = 1000;
    prm.ncycle = 51000;
    if SETUP >= 1 & SETUP <= 3
        prm.m = 3;
        if SETUP == 1
            prm.method = 'EnKF';
            prm.inflation = 1.35;
        elseif SETUP == 2
            prm.method = 'IEnKF';
            prm.inflation = 1.08;
        elseif SETUP == 3
            prm.method = 'IEKF';
            prm.inflation = 1.06;
        end
    elseif SETUP >= 4 & SETUP <= 8
        prm.m = 10;
        if SETUP == 4
            prm.method = 'EnKF';
            prm.inflation = 1.15;
        elseif SETUP == 5
            prm.method = 'IEnKF';
            prm.inflation = 1.02;
        elseif SETUP == 6
            prm.method = 'IEKF';
            prm.inflation = 1.06;
        elseif SETUP == 7
            prm.rotate = 1;
            prm.method = 'EnKF';
            prm.inflation = 1.04;
        elseif SETUP == 8
            prm.method = 'IEnKF';
            prm.inflation = 1.00;
         end
    end
elseif SETUP >= 9 & SETUP <= 16
    prm.cycle_nstep = 12;
    prm.obs_variance = 8;
    prm.nspinup = 1000;
    prm.ncycle = 101000;
    if SETUP >= 9 & SETUP <= 11
        prm.m = 3;
        if SETUP == 9
            prm.method = 'EnKF';
            prm.inflation = 1.08;
        elseif SETUP == 10
            prm.method = 'IEnKF';
            prm.inflation = 1.06;
        elseif SETUP == 11
            prm.method = 'IEKF';
            prm.inflation = 1.08;
        end
    elseif SETUP >= 12 & SETUP <= 16
        prm.m = 10;
        if SETUP == 12
            prm.method = 'EnKF';
            prm.inflation = 1.04;
        elseif SETUP == 13
            prm.method = 'IEnKF';
            prm.inflation = 1.02;
        elseif SETUP == 14
            prm.method = 'IEKF';
            prm.inflation = 1.08;
        elseif SETUP == 15
            prm.rotate = 1;
            prm.method = 'EnKF';
            prm.inflation = 1.04;
        elseif SETUP == 16
            prm.rotate = 1;
            prm.method = 'IEnKF';
            prm.inflation = 1.00;
        end
    end
end
