region='global';

figure
subplot(2,2,1)
set(gca,'fontsize',12)
load sampleerror_P40.mat
v=M;
plotRMSE1
text(1,14,'(a) without localization','fontsize',14)
ylabel('RMSE','fontsize',14)

subplot(2,2,2)
plotCR1
text(1,1.6,'(b) without localization','fontsize',14)
ylabel('CR','fontsize',14)

subplot(2,2,3)
load sampleerror_P40_localized.mat
plotlegend
plotRMSE1
text(1,14,'(c) with localization','fontsize',14)
ylabel('RMSE','fontsize',14)
xlabel('N','fontsize',14)

subplot(2,2,4)
plotCR1
text(1,1.6,'(d) with localization','fontsize',14)
ylabel('CR','fontsize',14)
xlabel('N','fontsize',14)

