
clear all
load Truth_N40_M80_P40.mat
M=40;
P=40;
obloc=1:P;
ob=ob(:,obloc);
R=R(obloc,obloc);
H=H(obloc,:);
d=1.2;
a=0.8;
F=[5 6 7 7.5 7.9 8];
lc=[0.01 0.1 1 10];
rc=0.0:0.1:1.0;
tau=100;
nk=1;

for i=1:length(F)
  [RMSE0_global(i,1,:),RMSE0_land(i,1,:),RMSE0_ocean(i,1,:),CR0(i,1,:,:),...
   lambda0(i,1,:,:),rcoef0(i,1,:,:),infx0(i,1,:,:,:)]=trial(a,1,d,F(i),dt,nk,nt,N,M,P,...
   u0,u1k,ob,obloc,R,H,false,'',0,0,'',1,0,false,0,tau);
  for j=1:length(lc)
    [RMSE1_global(i,j,:),RMSE1_land(i,j,:),RMSE1_ocean(i,j,:),CR1(i,j,:,:),...
     lambda1(i,j,:,:),rcoef1(i,j,:,:),infx1(i,j,:,:,:)]=trial(a,1,d,F(i),dt,nk,nt,N,M,P,...
     u0,u1k,ob,obloc,R,H,false,'anderson',ones(1,N),lc(j)*ones(1,N),'',1,0,false,0,tau);
  end
  for j=1:length(rc)
    [RMSE2_global(i,j,:),RMSE2_land(i,j,:),RMSE2_ocean(i,j,:),CR2(i,j,:,:),...
     lambda2(i,j,:,:),rcoef2(i,j,:,:),infx2(i,j,:,:,:)]=trial(a,1,d,F(i),dt,nk,nt,N,M,P,...
     u0,u1k,ob,obloc,R,H,false,'',0,0,'rtps',1,rc(j),false,0,tau);
  end
  [RMSE3_global(i,1,:),RMSE3_land(i,1,:),RMSE3_ocean(i,1,:),CR3(i,1,:,:),...
   lambda3(i,1,:,:),rcoef3(i,1,:,:),infx3(i,1,:,:,:)]=trial(a,1,d,F(i),dt,nk,nt,N,M,P,...
   u0,u1k,ob,obloc,R,H,false,'',0,0,'rtps_adaptive',1,0,false,0,tau);
end
save modelerror_adF_P40M40.mat

