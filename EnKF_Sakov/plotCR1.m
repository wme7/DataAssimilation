%figure
mean_cr0=nanmean(nanmean(CR0,3),4);
mean_cr1=nanmean(nanmean(CR1,3),4);
mean_cr2=nanmean(nanmean(CR2,3),4);
mean_cr3=nanmean(nanmean(CR3,3),4);

cmap=colormap('jet');
set(gca,'fontsize',10)
plot(v,ones(size(v)),'k:','linewidth',1);
hold on
plot(v,mean_cr0,'kx-','linewidth',1,'markersize',8);
for i=1:size(rc,2)
  plot(v,mean_cr2(:,i),'-','color',cmap(i*floor(size(cmap,1)/size(rc,2)),:),'linewidth',1,'markersize',8);
end
plot(v,mean_cr1(:,2),'k^','linewidth',1,'markersize',8);
plot(v,mean_cr1(:,3),'ks','linewidth',1,'markersize',8);
plot(v,mean_cr3,'ko','linewidth',1,'markersize',8);

%ylabel('consistency');
axis([min(v) max(v) 0 1.5]);
set(gca,'xtick',v,'xticklabel',M);
set(gca,'ytick',[0.1 0.5 1 1.2 1.5]);

