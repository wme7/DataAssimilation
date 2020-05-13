%load modelerror_F_P20M10_localized
%load sampleerror_P20_localized

set(gca,'fontsize',12)

plot(NaN,'k-');
hold on
plot(NaN,'^','color',[.5 .5 .5],'markersize',5);
plot(NaN,'s','color',[.5 .5 .5],'markersize',5);
plot(NaN,'x','color',[.5 .5 .5],'markersize',5);
plot(NaN,'ko','markersize',5);
legend('(color) RTPS','ACI 0.1','ACI 1','ACI 10','ACR','orientation','vertical')
colorbar('ytick',[1 33 65],'yticklabel',[0.0 0.5 1.0])

i=1;
plot(ones(1,N),'-','color',[.5 .5 .5])
cmap=colormap('jet');
for j=1:size(rc,2)
  plot(squeeze(nanmean(nanmean(infx2(i,j,:,:,:),3),4))','-','color',cmap(j*floor(size(cmap,1)/size(rc,2)),:))
end
plot(squeeze(nanmean(nanmean(infx1(i,2,:,:,:),3),4))','^','color',[.5 .5 .5],'markersize',5)
plot(squeeze(nanmean(nanmean(infx1(i,3,:,:,:),3),4))','s','color',[.5 .5 .5],'markersize',5)
plot(squeeze(nanmean(nanmean(infx1(i,4,:,:,:),3),4))','x','color',[.5 .5 .5],'markersize',5)
plot(squeeze(nanmean(infx3(i,1,1,:,:),4)),'ko','markersize',5)
hold off

axis([1 40 0.9 1.4])
xlabel('k')
ylabel('inflation')
set(gca,'ytick',[0.9 1 1.1 1.2 1.3])
