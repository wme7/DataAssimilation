function [x,y,z] = Lorenz63(rho, sigma, beta, X0, T, dt)
% Lorenz63 Function: solves explicitly the lorenz attractor of the prescribed
% values of parameters rho, sigma, beta
%
%       [x,y,z] = Lorenz63(rho, sigma, beta, X0, T, dt)
%
%       where: 
%
%       rho     - Rayleigh number
%       sigma   - Prandtl number
%       beta    - parameter
%       X0      - initial point (x0,y0,z0)
%       T       - time range [0,tEnd]
%       dt      - time step
%       [x,y,z] - output vectors of the strange attactor trajectories
%
% Example.
%        [x y z] = Lorenz63(28, 10, 8/3, [1,1,1], 1.0, 0.01); or
%        [x y z] = Lorenz63(28, 10, 8/3); or
%        [x y z] = Lorenz63();
%        plot3(x y z);
%
if nargin == 0
    rho  = 28;
    sigma= 10;    
    beta = 8/3;
    X0   = [0 1 1.05];
    T    = [0,25];
    dt   = 0.01;
elseif nargin<4
    X0   = [0 1 1.05];
    T    = [0,25];
    dt   = 0.01;
end

% Lorenz atractor function (typical values: rho=28; sigma=10; beta=8/3;)
%   x' = sigma*(y-x)
%   y' = x*(rho - z) - y
%   z' = x*y - beta*z
f = @(t,X) [sigma*(X(2)-X(1)), X(1)*(rho-X(3))-X(2), X(1)*X(2)-beta*X(3)];

t = (T(1):dt:T(end))';	% build time domain
dimX = length(X0);      % spatial size
dimT = length(t);       % temporal size
X = zeros(dimT,dimX);   % solution array 
X(1,:) = X0;            % Load initial condition

% Produce an explicit RK4 time integration
for n = 1:dimT-1
    k1 = f(t(n,:)       ,X(n,:)          );
    k2 = f(t(n,:)+0.5*dt,X(n,:)+0.5*dt*k1);
    k3 = f(t(n,:)+0.5*dt,X(n,:)+0.5*dt*k2);
    k4 = f(t(n,:)+  dt  ,X(n,:)+  dt*k3  );
    X(n+1,:)=X(n,:)+dt*(k1+2*k2+2*k3+k4)/6;
end

% Prepare function output
x=X(:,1); y=X(:,2); z=X(:,3);
end