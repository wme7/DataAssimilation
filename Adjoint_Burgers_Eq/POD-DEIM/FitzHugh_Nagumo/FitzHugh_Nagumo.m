close all
clear all
clf

ffun = @(v) (v.*(v-0.1).*(1-v)); %pre-define nonlinear function

global A C N eps dx E U Udeim Pdeim PP

%% Params
L     = 1;         % Parameters of the F-N system.
T     = 8;
eps   = 0.015;
b     = 0.5;
gamma = 2;
c     = 0.05;

N = 1024;          % spatial discretization full-order system
M = 1000;          % time discretization
ns= 100;           % # snapshots

%% Discretization in space and time
x=linspace(0,L,N);
dx =(L-0)/(N-1);
dt =(T-0)/(M-1); 

%% Solve Full Order System
v0 = 0*x';
w0 = 0*x';
y0 = [v0;w0]; %initial vector of unknowns

E = [eps*speye(N) zeros(N);
     zeros(N) speye(N)];
 
K = gallery('tridiag',N,-1,2,-1);
K(1,1)=1; K(end,end)=1; %Neumann BC

A = [-eps^2/dx^2*K -speye(N);
     b*speye(N) -gamma*speye(N)];
 
C = c*ones(2*N,1); 

%Use ode15s and define obtions
options = odeset('Mass',E,'MassSingular','no','Jacobian',@JacFN);

tic
[Tss,Yss] = ode15s(@rhs,linspace(0,T,ns),y0,options);  % compute snapshot matrix
[TOUT,YOUT] = ode15s(@rhs,linspace(0,T,M),y0,options); % compute accurate solution 
                                                       %  at equidistant time instances (M >> ns)
toc
Yss  = transpose(Yss);
YOUT = transpose(YOUT);

%% PLOT full order solution
[time, space] = meshgrid(TOUT,x);

figure(1)
Neval=1:N/32-1:N;
plot3(meshgrid(x(Neval),ones(M,1)),YOUT(Neval,:)',YOUT(N+Neval,:)','b:')
title('Phase-Space diagram')
xlabel('x') 
ylabel('v(x,t)') 
zlabel('w(x,t)')

%% POD-DEIM
[U1,S1,~] = svd(Yss(1:N,:)); 
[U2,S2,~] = svd(Yss(N+1:2*N,:)); 
F = ffun(Yss(1:N,:));
[U3,S3,~]=svd(F);

% The decay of singular values
figure(2)
semilogy(1:ns,diag(S1),'r+',1:ns,diag(S2),'+',1:ns,diag(S3),'k+');
ylim([1e-20 1e10]);
title('Decay of singular values')

% Solve reduced system
r=5;
U = blkdiag(U1(:,1:r),U2(:,1:r));
Udeim=U3(:,1:r);
[Pdeim,~]=deim(Udeim);              % DEIM algorithm
PP=Udeim*inv(Pdeim'*Udeim)*Pdeim';  % pre-compute projector

r1=size(U,2);
a0=zeros(r1,1);

tic
[~,aRed] = ode45(@rhsPODDEIM,linspace(0,T,M),a0);
toc
aRed=transpose(aRed);
YRed=U*aRed;

%% Plot reduced solution together with original system
figure(1)
hold on
plot3(meshgrid(x(Neval),ones(M,1)),YRed(Neval,:)',YRed(N+Neval,:)','r')

%% Numerical Analysis

PODerr= YOUT(1:N,:)-YRed(1:N,:);

figure(4)
surf(space,time,PODerr)                     % plot the error in space
shading interp

PODvec = reshape(transpose(PODerr),N*M,1);
L2errPOD      = dx*dt*norm(PODvec)          % error in L_2
MAXerrPOD     = norm(PODvec,inf)            % error in L_\infty