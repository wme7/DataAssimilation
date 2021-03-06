%figure
cmap=colormap('jet');
set(gca,'fontsize',10)
plot(ones(size(v)),'k:','linewidth',1);
hold on
plot(cr0,'kx-','linewidth',1,'markersize',8);
for i=1:size(rc,2)
  plot(cr2(:,i),'-','color',cmap(i*floor(size(cmap,1)/size(rc,2)),:),'linewidth',1,'markersize',8);
end
plot(cr1(:,2),'k^','linewidth',1,'markersize',8);
plot(cr1(:,5),'ks','linewidth',1,'markersize',8);
plot(cr3,'ko','linewidth',1,'markersize',8);

%ylabel('consistency');
axis([1 size(v,2) 0 1.5]);
set(gca,'xtick',1:size(v,2),'xticklabel',v);

%legend('1','no inflation',...
%'RTPS \alpha=0.1',...
%'RTPS \alpha=0.2',...
%'RTPS \alpha=0.3',...
%'RTPS \alpha=0.4',...
%'RTPS \alpha=0.5',...
%'RTPS \alpha=0.6',...
%'RTPS \alpha=0.7',...
%'RTPS \alpha=0.8',...
%'RTPS \alpha=0.9',...
%'ACI','ACR','location','NorthEastOutside')
