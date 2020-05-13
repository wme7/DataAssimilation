region='global';

figure
%subplot(2,2,1)
set(gca,'fontsize',12)
load sampleerror_P40.mat
v=[5 10 13 15 17 20 25 30 40]
plotRMSE1
text(5,14,'(a) without localization','fontsize',14)
ylabel('RMSE','fontsize',14)

%subplot(2,2,3)
%plotCR1
%text(5,1.6,'(b) without localization','fontsize',14)
%ylabel('CR','fontsize',14)
%xlabel('N','fontsize',14)
%
%subplot(2,2,2)
%load sampleerror_P40_localized.mat
%v=[4 5 6 7 8 9 10 12 14 16];
%plotRMSE1
%text(4,14,'(c) with localization','fontsize',14)
%ylabel('RMSE','fontsize',14)
%
%subplot(2,2,4)
%plotlegend
%plotCR1
%text(4,1.6,'(d) with localization','fontsize',14)
%ylabel('CR','fontsize',14)
%xlabel('N','fontsize',14)
%
