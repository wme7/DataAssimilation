
figure
subplot(3,1,1)
set(gca,'fontsize',12)
load sampleerror_P20_localized.mat
v=[4 5 6 7 8 9 10 12 14 16];
region='land';
plotlegend
plotRMSE
text(4,3.2,'(a) RMSE\_land','fontsize',14)
axis([min(v) max(v) 0 3])
set(gca,'ytick',[0.1 0.5 1 2 3]);

subplot(3,1,2)
set(gca,'fontsize',12)
region='ocean';
plotRMSE
text(4,3.2,'(b) RMSE\_ocean','fontsize',14)
axis([min(v) max(v) 0.9 3])
set(gca,'ytick',[1 2 3]);

subplot(3,1,3)
plotCR1
text(4,1.6,'(c) CR','fontsize',14)
axis([min(v) max(v) 0.5 1.5])
xlabel('N','fontsize',14)

