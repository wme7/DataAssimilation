prm.modelname = 'L40';
if 0 == 1
    prm.propagate = @L40_step;
else
    prm.propagate = @L40_step_c;
end
prm.generate = @L40_generate;
prm.n = 40;
prm.m = 25;
plotfunc = @semilogy;
prm.collapseratio = 0.2;

% set the experiment here

SETUP = 2;

% set the output to your convenience

prm.doplot = 0;
prm.verbose = 0;
prm.calc_norms = 1;

% do not edit below

if SETUP >= 1 & SETUP <= 8
    prm.observe = @L40_observe_1;
    prm.p = length(prm.observe([1 : prm.n], 0));
    prm.cycle_nstep = 12;
    prm.obs_variance = 1;
    prm.ncycle = 51000;
    prm.nspinup = 1000;
    if SETUP >= 1 & SETUP <= 3
        prm.m = 25;
        if SETUP == 1
            prm.method = 'EnKF';
            prm.inflation = 1.80;
        elseif SETUP == 2
            prm.method = 'IEnKF';
            prm.inflation = 1.20;
        elseif SETUP == 3
            prm.method = 'IEKF';
            prm.inflation = 1.50;
        end
    elseif SETUP >= 4 & SETUP <= 8
        prm.m = 60;
        if SETUP == 4
            prm.method = 'EnKF';
            prm.inflation = 1.25;
        elseif SETUP == 5
            prm.method = 'IEnKF';
            prm.inflation = 1.15;
        elseif SETUP == 6
            prm.method = 'IEKF';
            prm.inflation = 1.50;
        elseif SETUP == 7
            prm.method = 'EnKF';
            prm.rotate = 1;
            prm.inflation = 1.20;
        elseif SETUP == 8
            prm.method = 'IEnKF';
            prm.rotate = 1;
            prm.inflation = 1.10;
        end
    end
elseif SETUP >= 9 & SETUP <= 14
    prm.observe = @L40_observe_5;
    prm.p = length(prm.observe([1 : prm.n], 0));
    prm.obs_variance = 0.001;
    prm.ncycle = 12625;
    prm.nspinup = 125;
    if SETUP >= 9 & SETUP <= 11
        prm.cycle_nstep = 7;
        if SETUP == 9
            prm.method = 'EnKF';
            prm.inflation = 1.08;
        elseif SETUP == 10
            prm.method = 'IEnKF';
            prm.inflation = 1.02;
        elseif SETUP == 11
            prm.method = 'IEKF';
            prm.inflation = 1.10;
        end
    elseif SETUP >= 12 & SETUP <= 14
        prm.cycle_nstep = 6;
        if SETUP == 12
            prm.method = 'EnKF';
            prm.inflation = 1.04;
        elseif SETUP == 13
            prm.method = 'IEnKF';
            prm.inflation = 1.02;
        elseif SETUP == 14
            prm.method = 'IEKF';
            prm.inflation = 1.06;
        end
    end
elseif SETUP >= 15 & SETUP <= 23
    prm.observe = @L40_observe_4;
    prm.p = length(prm.observe([1 : prm.n], 0));
    prm.obs_variance = 0.001;
    prm.ncycle = 12625;
    prm.nspinup = 125;
    if SETUP >= 15 & SETUP <= 17
        prm.cycle_nstep = 12;
        if SETUP == 15
            prm.method = 'EnKF';
            prm.inflation = NaN;
        elseif SETUP == 16
            prm.method = 'IEnKF';
            prm.inflation = 1.06;
        elseif SETUP == 17
            prm.method = 'IEKF';
            prm.inflation = NaN;
        end
    elseif SETUP >= 18 & SETUP <= 20
        prm.cycle_nstep = 10;
        if SETUP == 18
            prm.method = 'EnKF';
            prm.inflation = NaN; % 1.15
        elseif SETUP == 19
            prm.method = 'IEnKF';
            prm.inflation = 1.02;
        elseif SETUP == 20
            prm.method = 'IEKF';
            prm.inflation = 1.15;
        end
    elseif SETUP >= 21 & SETUP <= 23
        prm.cycle_nstep = 8;
        if SETUP == 21
            prm.method = 'EnKF';
            prm.inflation = 1.10;
        elseif SETUP == 22
            prm.method = 'IEnKF';
            prm.inflation = 1.02;
        elseif SETUP == 23
            prm.method = 'IEKF';
            prm.inflation = 1.08;
        end
    end
elseif SETUP >= 24 & SETUP <= 31
    prm.observe = @L40_observe_1;
    prm.p = length(prm.observe([1 : prm.n], 0));
    prm.obs_variance = 20;
    prm.cycle_nstep = 1;
    prm.calc_norms = 0;
    prm.ncycle = 101000;
    prm.nspinup = 1000;
    if SETUP >= 24 & SETUP <= 26
        prm.m = 25;
        if SETUP == 24
            prm.method = 'EnKF';
            prm.inflation = 1.04;
        elseif SETUP == 25
            prm.method = 'IEnKF';
            prm.inflation = 1.04;
        elseif SETUP == 26
            prm.method = 'IEKF';
            prm.inflation = 1.10;
        end
    elseif SETUP >= 27 & SETUP <= 31
        prm.m = 60;
        if SETUP == 27
            prm.method = 'EnKF';
            prm.inflation = 1.02;
        elseif SETUP == 28
            prm.method = 'IEnKF';
            prm.inflation = 1.02;
        elseif SETUP == 29
            prm.method = 'IEKF';
            prm.inflation = 1.10;
        elseif SETUP == 30
            prm.method = 'EnKF';
            prm.rotate = 1;
            prm.inflation = 1.03;
        elseif SETUP == 31
            prm.method = 'IEnKF';
            prm.rotate = 1;
            prm.inflation = 1.02;
        end
    end
end
if 1 == 0
    prm.nspinup = 200;
    prm.ncycle = 200 + prm.nspinup;
end
