function [] = stats_print(prm, stats)
    
    if ~exist('nbad', 'var')
        nbad = 0;
    end
    if isfield(prm, 'report_custom')
        prm.report_custom(stats);
    end
    good = find(isfinite(stats.rmse_a));
    fprintf('  mean analysis rmse = %.3g\n', mean(stats.rmse_a(good)));
    fprintf('  mean analysis spread = %.3g\n', mean(stats.spread_a(good)));
    fprintf('  mean forecast rmse = %.3g\n', mean(stats.rmse_f(good)));
    fprintf('  mean forecast spread = %.3g\n', mean(stats.spread_f(good)));
    fprintf('  mean # iterations = %.3g\n', mean(stats.niter(good)));
    if isfield(stats, 'normT')
        fprintf('  mean ||T^-1|| * ||A_1^f|| / ||A_2^f|| = %.4g\n', mean(stats.normT ./ (stats.normA2_f ./ stats.normA1_f)));
    end

    return
    