function [stats] = stats_merge(stats, stats_now)
    
    if isempty(stats)
        stats = stats_now;
        return
    end
    
    stats.rmse_f = [stats.rmse_f; stats_now.rmse_f];
    stats.rmse_a = [stats.rmse_a; stats_now.rmse_a];
    stats.rmse = [stats.rmse; stats_now.rmse];
    stats.spread_f = [stats.spread_f; stats_now.spread_f];
    stats.spread_a = [stats.spread_a; stats_now.spread_a];
    stats.kurt_f = [stats.kurt_f; stats_now.kurt_f];
    stats.kurt_a = [stats.kurt_a; stats_now.kurt_a];
    stats.inc = [stats.inc; stats_now.inc];
    stats.dy = [stats.dy; stats_now.dy];
    stats.niter = [stats.niter; stats_now.niter];
    if isfield(stats, 'normA1_f')
        stats.normA1_f = [stats.normA1_f; stats_now.normA1_f];
        stats.normA2_f = [stats.normA2_f; stats_now.normA2_f];
        stats.normT = [stats.normT; stats_now.normT];
        stats.normHA_f = [stats.normHA_f; stats_now.normHA_f];
        stats.normHA_a = [stats.normHA_a; stats_now.normHA_a];
    end
    if isfield(stats, 'nonl')
        stats.nonl = [stats.nonl; stats_now.nonl];
    end

    return
    