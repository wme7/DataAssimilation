%load M.mat
%v=M;
%figure
m_rmse0=nanmedian(eval(['RMSE0_' region ]),3);
m_rmse1=nanmedian(eval(['RMSE1_' region ]),3);
m_rmse2=nanmedian(eval(['RMSE2_' region ]),3);
m_rmse3=nanmedian(eval(['RMSE3_' region ]),3);

cmap=colormap('jet');
set(gca,'fontsize',10)
semilogy(v,sqrt(mean(diag(R)))*ones(size(v)),'-','color',[.5 .5 .5],'linewidth',1,'markersize',8);
hold on
semilogy(v,m_rmse0,'kx-','linewidth',1,'markersize',8);
for i=1:size(rc,2)
  semilogy(v,m_rmse2(:,i),'-','color',cmap(i*floor(size(cmap,1)/size(rc,2)),:),'linewidth',1,'markersize',8);
end
semilogy(v,m_rmse1(:,2),'k^','linewidth',1,'markersize',8);
semilogy(v,m_rmse1(:,3),'ks','linewidth',1,'markersize',8);
semilogy(v,m_rmse3,'ko','linewidth',1,'markersize',8);

%ylabel('m_rmse');
axis([min(v) max(v) 0.1 10]);
set(gca,'xtick',v,'xticklabel',M);

