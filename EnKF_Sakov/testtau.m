
clear all
load Truth_N40_M80_P40.mat
M=[20 40];
P=40;
obloc=1:P;
ob=ob(:,obloc);
R=R(obloc,obloc);
H=H(obloc,:);
d=1;
a=1;
F=8;
tau=[1 100]; %[1 10 100 1000];
nk=1;

for i=1:length(M)
for j=1:length(tau)
  [RMSE3_global(i,j,:),RMSE3_land(i,j,:),RMSE3_ocean(i,j,:),CR3(i,j,:,:),...
   lambda3(i,j,:,:),rcoef3(i,j,:,:),infx3(i,j,:,:,:)]=trial(a,1,d,F,dt,nk,nt,N,M(i),P,...
   u0,u1k,ob,obloc,R,H,false,'',0,0,'rtps_adaptive',1,0,false,0,tau(j));
save 2.mat
end
end

