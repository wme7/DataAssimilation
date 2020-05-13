function [rmse_global,rmse_land,rmse_ocean,cr,lambda,r_coef,infx]=trial(a,b,d,F,dt,nk,nt,N,M,P,...
                    u0,u1k,ob,obloc,R,H,...
                    inf_additive,inflation_type,inf,inf_var,...
                    relaxation_type,lambda0,r_coef0,...
                    localize,l_coef,tau)
Q=0.01*eye(N);
r_coef=nan(1,1,nk,nt+2);
lambda=nan(1,1,nk,nt+2);
r_coef(1,1,:,1)=r_coef0;
lambda(1,1,:,1)=lambda0;
rmse=nan(1,1,nk);
cr=nan(1,1,nk,nt);
infx=nan(1,1,nk,nt,N);

for k=1:nk
%  figure('position',[0 0 800 200])

  u1=nan(nt+1,M,N);
  u2=nan(nt+1,M,N);
  for m=1:M
    u1(1,m,:)=u1k(k,m,:);
  end

  for t=1:nt
    if(max(u2(:))>1e2)
      display(['model blow up after ' num2str(t) ' steps'])
      break
    end
    if(t>0)
      if(inf_additive & t>1)
        [Q,R,status]=inflation_additive(N,M,P,H,u1(t-1,:,:),u1(t,:,:),u2(t-1,:,:),ob(t-1,:),ob(t,:),Q,R);
        if(status>0)
          break
        end
      end

      [u2(t,:,:),lambda(1,1,k,t+1),r_coef(1,1,k,t+1),cr(1,1,k,t),infx(1,1,k,t,:)]=enkf(u1(t,:,:),...
                   N,P,M,H,ob(t,:),obloc,R,...
                   Q,inf_additive,inflation_type,inf,inf_var,...
                   relaxation_type,lambda(1,1,k,t),r_coef(1,1,k,t),...
                   localize,l_coef,tau);
    else
      u2(t,:,:)=u1(t,:,:);
    end
    u1(t+1,:,:)=l96(u2(t,:,:),dt,N,a,b,d,F);

%  set(gca,'fontsize',12)
%  plot(squeeze(u1(t,:,:))','c-','linewidth',1);
%  hold on
%  plot(squeeze(mean(u1(t,:,:),2)),'b-','linewidth',2);
%  plot(u0(t,:),'k-','linewidth',2);
%  plot(squeeze(u2(t,:,:))','y-','linewidth',1);
%  plot(squeeze(mean(u2(t,:,:),2)),'r-','linewidth',2);
%  plot(obloc,ob(t,:),'k*')
%  plot(squeeze(infx(1,1,k,t,:)))
%  hold off
%  axis([1 N -10 10])
%  xlabel('k')
%  pause(0.01)
%if(t>1)
%plot(squeeze((std(u1(t,:,:),0,2)-std(u2(t,:,:),0,2))./(std(u2(t,:,:),0,2))),'k')
%plot(squeeze(std(u1(t,:,:),0,2)./std(u2(t-1,:,:),0,2)),'k')
%axis([1 N 0 2])
%pause(0.1)
%end

  end

  udiff1=mean(u2(1:nt,:,:),2)-u0(1:nt,1,:);
  udiff2=mean(u2(1:nt,:,obloc),2)-u0(1:nt,1,obloc);
  udiff3=mean(u2(1:nt,:,21:40),2)-u0(1:nt,1,21:40);
  rmse_global(1,1,k)=sqrt(sum(sum(udiff1(:).*udiff1(:)))/nt/N);
  rmse_land(1,1,k)=sqrt(sum(sum(udiff2(:).*udiff2(:)))/nt/N);
  rmse_ocean(1,1,k)=sqrt(sum(sum(udiff3(:).*udiff3(:)))/nt/N);

end

