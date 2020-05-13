clear;
SEMILOG = 1;

METHODS = {'IEnKF', 'IEKF', 'EnKF'};
INFLS = [0.08 0.10 0.10 0.10 0.12 0.12 0.12 0.12 0.14 0.14 0.14 0.18 0.20 0.20 0.24 0.26;
         0.34 0.42 0.40 0.40 0.40 0.40 0.40 0.40 0.40 0.42 0.44 0.52 0.62 0.72 0.72 0.72;
         0.08 0.14 0.20 0.26 0.32 0.34 0.34 0.34 0.40 0.48 0.48 0.50 0.50 0.52 0.52 0.52];

[NM, NT] = size(INFLS);
rmse = zeros(size(INFLS));
rmse(:, :) = NaN;

for m = 1 : NM
    method = char(METHODS{m});
    for i = 1 : NT
        infl = INFLS(m, i);
        cycle_nstep = i;
        fname = sprintf('out/%s-%d-%.2f.mat', method, cycle_nstep, infl);
        if exist(fname, 'file')
            fprintf('  %s +\n', fname);
            load(fname);
            rmse(m, i) = mean(stats.rmse_a);
        else
            fprintf('  %s -\n', fname);
        end
    end
end

figure(1);
if ~SEMILOG
    plot(rmse','.-');
else
    semilogy(rmse','.-');
end
xlim([1 NT]);
ylim([0 1.2]);
xlabel('T');
ylabel('RMSE');
legend(METHODS, 'location', 'northwest');
grid;
line(xlim, [1 1], 'linestyle', '--', 'color', 'k');

if SEMILOG
    set(gca, 'YTickLabelMode', 'manual');
    set(gca, 'YTick', [0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2]);
    set(gca, 'YTickLabel', ['0.2'; '0.3'; '0.4'; '0.5'; '0.6'; '0.7'; '0.8'; '0.9'; '1.0'; '   '; '1.2']);
end
