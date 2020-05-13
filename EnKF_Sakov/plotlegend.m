set(gca,'fontsize',10)
plot(NaN,'-','color',[.5 .5 .5],'linewidth',1,'markersize',8);
hold on
plot(NaN,'kx-','linewidth',1,'markersize',8);
plot(NaN,'k-','linewidth',1,'markersize',8);
plot(NaN,'k^','linewidth',1,'markersize',8);
plot(NaN,'ks','linewidth',1,'markersize',8);
plot(NaN,'ko','linewidth',1,'markersize',8);
legend('obs noise','no inflation','(color) RTPS','ACI 0.1','ACI 1','ACR','orientation','horizontal')
colorbar('ytick',[1 33 65],'yticklabel',[0.0 0.5 1.0])
