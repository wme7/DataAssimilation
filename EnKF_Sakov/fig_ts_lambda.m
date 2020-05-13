load tau_F_P20M10_localized

figure
subplot(2,1,1)
set(gca,'fontsize',12)
i=4;
plot(1:200,squeeze(lambda3(i,1,1:200)),'b'); hold on
plot(250:450,squeeze(lambda3(i,1,4800:5000)),'b');
plot(1:200,squeeze(lambda3(i,2,1:200)),'r');
plot(250:450,squeeze(lambda3(i,2,4800:5000)),'r');
plot(1:200,squeeze(lambda3(i,3,1:200)),'k'); 
plot(250:450,squeeze(lambda3(i,3,4800:5000)),'k'); hold off
axis([0 450 0 2])
text(1,2.2,'(a) perfect model','fontsize',14)
text(220,1.0,'...','fontsize',15)
set(gca,'xticklabel',[0:50:200 4800:50:5000])
ylabel('\lambda','fontsize',14)

subplot(2,1,2)
set(gca,'fontsize',12)
i=1;
plot(1:200,squeeze(lambda3(i,1,1:200)),'b'); hold on
plot(250:450,squeeze(lambda3(i,1,4800:5000)),'b');
plot(1:200,squeeze(lambda3(i,2,1:200)),'r');
plot(250:450,squeeze(lambda3(i,2,4800:5000)),'r');
plot(1:200,squeeze(lambda3(i,3,1:200)),'k'); 
plot(250:450,squeeze(lambda3(i,3,4800:5000)),'k'); 
axis([0 450 0 2])
text(1,2.2,'(b) model error F=5','fontsize',14)
text(220,1.0,'...','fontsize',15)
xlabel('time step','fontsize',14)
ylabel('\lambda','fontsize',14)
set(gca,'xticklabel',[0:50:200 4800:50:5000])
plot([150 180],[0.2 0.2],'b','linewidth',2);
plot([250 280],[0.2 0.2],'r','linewidth',2);
plot([350 380],[0.2 0.2],'k','linewidth',2);
text(190,0.2,'\tau=1','fontsize',14);
text(290,0.2,'\tau=10','fontsize',14);
text(390,0.2,'\tau=100','fontsize',14);


