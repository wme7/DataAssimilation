function [x_true, E] = L40_generate(prm)

    if prm.obs_variance >= 1
        load('L40_spinned_1.mat');
    else
        load('L40_spinned_0.1.mat');
    end
    mm = size(E, 2);
    if prm.m > mm
        error('\n  error: L40_generate(): m (= %d) > m_max (= %d)\n', prm.m, mm);
    end
    rand_state = rand('state');
    randn_state = randn('state');
    pool = shuffle(mm);
    E = E(:, pool(1 : prm.m));
    rand('state', rand_state);
    randn('state', randn_state);
    
    randn(mm, 1); % to avoid correlation between observations and ensemble
    
    return
