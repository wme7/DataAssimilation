close all
clear all
clf 

global nu tol_newton max_newton h dt M B C Nt N k1 k2 ...
    time space plotOn Upod Bred Cred Fred

addpath('./FullBurgers');
flag = 2;

if flag == 1
    addpath('./POD');
    disp('POD')
else
    addpath('./PODDEIM');
    disp('POD-DEIM')
end


%% Params
N  = 80;                % # grid point in space
Nt = 80;                % # grid point in time
L  = 1;                 % x \in [0,L]
T  = 1;                 % t \in [0,T]

x = linspace(0,L,N+1);           % Spacial discretization
t = linspace(0,T,Nt);            % time discretization
[time, space] = meshgrid(t',x);  % for 3D plotting
h = L/N;                         % constant step size
dt= T/(Nt-1);                    % constant time step

nu    = 0.01;           % viscosity parameter

tol_newton = 10e-3;     % for Newton's method
max_newton = 20;        % max. number of inner Newton iterations

k1 = 9;                 % # POD basis of y
k2 = 13;                % # DEIM basis

plotOn = 1;             % plotting and output on (1) or off(0)


%% Initialize the state y_0
y0 = [ones(floor((N-1)/2),1);zeros(ceil((N-1)/2),1)];

% Pre-compute stiffness and mass matrix for full-order Burgers eqn
M = (h/6)*gallery('tridiag',ones(N-2,1),4*ones(N-1,1),ones(N-2,1));
B = gallery('tridiag',-0.5*ones(N-2,1),zeros(N-1,1),0.5*ones(N-2,1));
C = (1/h)*gallery('tridiag',-ones(N-2,1),2*ones(N-1,1),-ones(N-2,1));


%% Solve Full-order system
tStart = tic;
Yfull = Burgers(y0);
tElapsed = toc(tStart)
Yfinal=Yfull;


%% Setup POD
if flag == 1
    tStart = tic;
    R = chol(M);
    [Uy,Sy,~]=svd(R*Yfull);
    Upod=R\Uy(:,1:k1);

    if plotOn ==1
       figure(55)
       index=1:min(N-1,Nt);
       semilogy(index,diag(Sy(index,index)),'x'); %investigate basis dimensions
    end


    % Pre-compute matrices
    Cred  = Upod'*C*Upod;
    Bred  = Upod'*B;

    % obtain reduced variables
    y0_red = Upod'*M*y0;
    tElapsed_setup = toc(tStart)  

else

    %% Setup POD-DEIM
    tStart = tic;
    % SVD of snapshot matrices + DEIM
    R = chol(M);
    [Uy,Sy,~]=svd(R*Yfull);
    Upod=R\Uy(:,1:k1);
    
    [Uf,Sf,~]=svd(Yfull.^2);
    Udeim=Uf(:,1:k2);
    [~,ind]=deim(Udeim);

    if plotOn ==1
       figure(55)
       index=1:min(N-1,Nt);
       semilogy(index,diag(Sy(index,index)),'x',index,diag(Sf(index,index)),'rx'); %investigate basis dimensions
    end
    
    % Pre-compute matrices
    Bred  = Upod'*B*Udeim/Udeim(ind,:);
    Cred  = Upod'*C*Upod;
    Fred  = Upod(ind,:);
    
    % obtain reduced variables
    y0_red = Upod'*M*y0;
    tElapsed_setup = toc(tStart) 
end
    
%% Solve red system
tStart = tic;
Y_red = Burgers_Reduced(y0_red);
tElapsed_red = toc(tStart)   
    

Yapprox = Upod*Y_red;

%% Plot full-order system and reduced model in one figure

figure(5)
subplot(1,2,1);
mesh(time, space, [zeros(1,Nt);Yfull;zeros(1,Nt)]); %add D-BC to Yfull
xlabel('t')
ylabel('x')
zlabel('y(t,x)'); zlim([0 1.5])
title(['Full Order System, n = ' num2str(N)])
subplot(1,2,2);
mesh(time, space, [zeros(1,Nt);Yapprox;zeros(1,Nt)]); %add D-BC to Yapprox
xlabel('t')
ylabel('x')
zlabel('y_k(t,x)'); zlim([0 1.5])
title(['# POD basis = ' num2str(k1) '     # DEIM basis = ' num2str(k2)])

%% Some numerical analysis

SpeedUp = tElapsed/(tElapsed_red+tElapsed_setup)
SpeedUp_solver = tElapsed/(tElapsed_red)

error = L2norm(Yfinal,Yapprox)


err   = Yfinal-Yapprox;
errL2 = sqrt(h*dt*sum(sum(err.^2)));
relErrL2= errL2/sqrt(h*dt*sum(sum(Yfinal.^2)))


%% remove paths for next run
warning off
rmpath('./FullBurgers');
rmpath('./POD');
rmpath('./PODDEIM');

