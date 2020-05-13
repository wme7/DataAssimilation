function [x_true, E] = L3_generate(prm)

    load('L3_samples', 'S');
    x_true = [1.508870; 21.531271; 25.46091];
    rand_state = rand('state');
    randn_state = randn('state');
    S = S(:, 1 : 1000);
    pool = shuffle(1000);
    E = S(:, pool(1 : prm.m));
    rand('state', rand_state);
    randn('state', randn_state);
    
    randn(100, 1); % to avoid correlation between observations and ensemble
                   % (assuming m <= 100)

    return
