%%
% Main file solving the optimal control problem
%
% \min_u \frac{1}{2} \int_0^T \int_0^L [y(t,x) - z(t,x)]^2 + \omega u^2(t,x) dx dt,
%
% where y is a solution to the nonlinear Burgers' equation
%
% y_t + ( \frac{1}{2}y^2 - \nu y_x )_x = f + u,  (x,t) \in (0,L) \times (0,T), 
% y(t,0) = y(t,L) = 0,  t \in (0,T),
% y(0,x) = y_0(x),  x \in (0,L).
%
% The optimization is done on a full-order, POD-, and POD-DEIM-reduced
% model using the optimization algorithms Newton-type, BFGS and SPG.
%
%%

close all
clear all
clf 

global nu tol_newton max_newton h dt omega M B C Nt N k1 k2...
    time space plotOn Upod iter_tmp

addpath('./FullBurgers');
%addpath('./POD');
addpath('./PODDEIM');

flag = 3;

if flag == 1 
    addpath('./SPG');
elseif flag == 2
    addpath('./BFGS','./BFGS/immoptibox/');
elseif flag == 3
    addpath('./Newton-type');
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
omega = 0.005;          % control penalty

tol_newton = 10e-3;     % for Newton's method
max_newton = 20;        % max. number of inner Newton iterations

k1 = 9;                 % # POD basis of y
k2 = 13;                % # DEIM basis

plotOn = 1;             % plotting and output on (1) or off(0)

%% desired state z, does not depend on time
z  = [ones(floor((N-1)/2),1);zeros(ceil((N-1)/2),1)];
Z  = kron(ones(1,Nt),z);

%% Initialize the state y_0
y0 = [ones(floor((N-1)/2),1);zeros(ceil((N-1)/2),1)];

% Pre-compute stiffness and mass matrix for full-order Burgers eqn
M = (h/6)*gallery('tridiag',ones(N-2,1),4*ones(N-1,1),ones(N-2,1));
B = gallery('tridiag',-0.5*ones(N-2,1),zeros(N-1,1),0.5*ones(N-2,1));
C = (1/h)*gallery('tridiag',-ones(N-2,1),2*ones(N-1,1),-ones(N-2,1));

%% Initialize the control u_0 = 0
u0 = zeros((N-1)*Nt,1);

%% Optimization - full order

tStart = tic;
if flag == 1 % ... using SPG method
    disp('SPG - full Model')
    
    Opts = struct('x0',u0,'n',(N-1)*Nt);
    xl = -2;
    xr = 2;
    Opts.Proj = @(x)min(max(x,xl),xr); 
    
    [u,INFO] = spg_cons(@(u)SPGfun(u, z, y0), Opts);
end

if flag == 2 % ... using BFGS method
    Opts = [1 10e-8 10e-6 40];
    iter_tmp = 0;
   
    u = ucminf(@(u)BFGSfun(u, z, y0), u0,Opts);
end


if flag == 3 % ...using Newton-eqn + adjoints
    Opts = struct('max_opt',30,'tol_opt_f',10e-8,'tol_opt_gradf',10e-9,'tol_line',10e-4,'min_alpha',10e-8,'eta',10e-2);

    [u,~] = minNewtonCG(u0,y0,z,Opts);
end
tElapsed = toc(tStart)

Yfinal = Burgers(reshape(u,N-1,Nt),y0);
costFull = costFunL2( u, Yfinal, z )


%% Optimization - reduced model
tStart = tic;

U = reshape(u0,N-1,Nt);
Yapprox = Burgers(U,y0);

red_iter = 0;
while (L2norm(Yapprox,Z) > 0.5 && red_iter < 5)
    
    [z_red, y0_red] = setUp_redModel(Yapprox, z, y0);

    if flag == 1 % ... using SPG
        disp('SPG - reduced Model')
        
        Opts = struct('x0',u0,'n',(N-1)*Nt);
        xl = -2;
        xr = 2;
        Opts.Proj = @(x) min(max(x,xl),xr); 

        u = spg_cons(@(u)SPGfun_Reduced(u, z_red, y0_red), Opts);
    end
       
    if flag == 2 % ... using BFGS
        iter_tmp = 0;
        
        u = ucminf(@(u)BFGSfun_Reduced(u, z_red, y0_red), u0,Opts);    
    end

    if flag == 3 % ... using Newton-eqn + adjoints
        
        Opts = struct('max_opt',7,'tol_opt_f',10e-8,'tol_opt_gradf',10e-9,'tol_line',10e-4,'min_alpha',10e-8,'eta',10e-2);
        u = minNewtonCG_Reduced(u,y0_red,z_red,Opts); 
    end

    Y_red = Burgers_Reduced(reshape(u,N-1,Nt),y0_red);
    Yapprox = Upod*Y_red;
   
    red_iter = red_iter + 1;
end
tElapsed_all = toc(tStart)


%% Some numerical analysis

SpeedUp = tElapsed/tElapsed_all
costFullfromRed = costFunL2( u, Yapprox, z )
costFull

err   = Yfinal-Yapprox;
errL2 = sqrt(h*dt*sum(sum(err.^2)));
relErrL2= errL2/sqrt(h*dt*sum(sum(Yfinal.^2)))

%% remove paths for next run
warning off
rmpath('./FullBurgers');
rmpath('./POD');
rmpath('./PODDEIM');
rmpath('./SPG');
rmpath('./BFGS','./BFGS/immoptibox/');
rmpath('./Newton-type');