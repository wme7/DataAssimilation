set(gca,'DataAspectRatio',[25 5 1]);
xlim([100 201]);
set(gca, 'xtick', 100 : 10 : 200);
set(gca, 'xticklabel', [0 : 10 : 100]);
grid;
set(gca, 'ygrid', 'off');
title('');
[a,b,c,d] = legend('', 'forecast', 'last iteration','','', '', '', '', '1st iteration', '2nd iteration', 'location', 'southwest');
legend(c([2 9 10 3]), d([2 9 10 3]), 'location', 'southwest');
set(gcf, 'papersize', [6 3], 'paperunits', 'inches');
print -depsc2 ienkf-100steps.ps;