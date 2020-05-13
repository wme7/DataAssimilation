clear all

aTrue=1;
b=1;
dTrue=1;
FTrue=8;
N=40;
M=80;
P=40;
dt=0.05;
nk=1;
nt=100000;
obloc=1:P;
%obloc=obs_location(N,P,'regular1');
H=h_operator(N,P,obloc);
R=1.0*eye(P);

u0=nan(nt+1,1,N);
u0(1,1,:)=randn(1,1,N);
for t=1:nt
  u0(1,:,:)=l96(u0(1,:,:),dt,N,aTrue,b,dTrue,FTrue);
end
for k=1:nk
for m=1:M
  u1k(k,m,:)=u0(1,1,:)+0.1*randn(1,1,N);
end
end
for t=1:nt+1
  ob(t,:)=H*squeeze(u0(t,1,:))+sqrt(R)*randn(P,1);
  u0(t+1,:,:)=l96(u0(t,:,:),dt,N,aTrue,b,dTrue,FTrue);
end
save Truth_N40_M80_P40.mat
