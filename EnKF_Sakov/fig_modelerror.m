region='global';

figure
subplot(3,2,1)
set(gca,'fontsize',12)
load save2/modelerror_F_P40M40.mat
v=F;
plotRMSE1
text(1,14,'(a) a=1; d=1','fontsize',14)
ylabel('RMSE','fontsize',14)

%subplot(3,2,2)
%plotCR1
%text(1,1.6,'(b) a=1; d=1','fontsize',14)
%ylabel('CR','fontsize',14)
%
%subplot(3,2,3)
%load modelerror_dF_P40M40.mat
%plotRMSE1
%text(1,14,'(c) a=1; d=1.2','fontsize',14)
%ylabel('RMSE','fontsize',14)
%
%subplot(3,2,4)
%plotCR1
%text(1,1.6,'(d) a=1; d=1.2','fontsize',14)
%ylabel('CR','fontsize',14)
%
%subplot(3,2,5)
%load modelerror_aF_P40M40.mat
%plotlegend
%plotRMSE1
%text(1,14,'(e) a=0.8; d=1','fontsize',14)
%ylabel('RMSE','fontsize',14)
%xlabel('F','fontsize',14)
%
%subplot(3,2,6)
%plotCR1
%text(1,1.6,'(f) a=0.8; d=1','fontsize',14)
%ylabel('CR','fontsize',14)
%xlabel('F','fontsize',14)
