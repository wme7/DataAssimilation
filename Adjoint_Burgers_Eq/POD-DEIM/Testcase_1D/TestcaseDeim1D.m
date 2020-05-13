close all
clear all
clf

%% Params
mu=3.1;       % parameter of interest
n=100;        % # spatial discretization
ns=51;        % # snapshots for \mu
m=10;         % # DEIM dimension

%% Define the nonlinear function f
x=linspace(-1,1,n);
f=(1-x).*cos(3*pi*mu.*(x+1)).*exp(-mu*(1+x));

%% DEIM
muVec=linspace(1,pi,ns);
F=zeros(n,ns);
for ii=1:ns %obtain snapshot matrix
    F(:,ii)= (1-x).*cos(3*pi*muVec(ii).*(x+1)).*exp(-muVec(ii)*(1+x));
end
[U,S,V]=svd(F);
U=U(:,1:m);
[P,IND] = deim(U);
fred=U*inv((transpose(P)*U))*transpose(P)*f';

%% Plot - Decay of singular values
figure(1)
semilogy(diag(S),'x')
title(['Sing values of ' num2str(ns) ' snapshots'])

%% Plot - POD bases and location of DEIM points
figure(2)
plot(x,U(:,1),x,U(:,2),'g',x,U(:,3),'r',x,U(:,4),'m',x,U(:,5),'k',x,U(:,6),'c')
hold on
plot(x(IND(1:m)),0,'bo','MarkerFaceColor','k')
title(['POD bases and DEIM indices (1-' num2str(m) ')'])

%% Plot - Full order and reduced model
figure(3)
plot(x,f,x,fred,'rx')
ylim([-2 2])
title(['\mu = ' num2str(mu)])