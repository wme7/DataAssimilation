function [x,y,z] = Lorenz84(a,b,F,G,X0,T,dt)
% Lorenz84 Function: solves explicitly the lorenz attractor of the prescribed
% values of parameters rho, sigma, beta
%
%       [x,y,z] = Lorenz84(rho, sigma, beta, X0, T, dt)
%
%       where: 
%
%       a,b,F,G - parameters
%       X0      - initial point (x0,y0,z0)
%       T       - time range [0,tEnd]
%       dt      - time step
%       [x,y,z] - output vectors of the strange attactor trajectories
%
% Example.
%        [x y z] = Lorenz84(0.25, 4, 8.0, 1.0, [1,1,1], 1.0, 0.01); or
%        [x y z] = Lorenz84(0.25, 4, 8.0, 1.0); or
%        [x y z] = Lorenz84();
%        plot3(x y z);
%
if nargin == 0
    a  = 0.25;
    b  = 4.0;    
    F  = 8.0;
    G  = 1.0;
    X0 = [0 0 0];
    T  = [0,25];
    dt = 0.01;
elseif nargin<4
    X0 = [0 0 0];
    T  = [0,25];
    dt = 0.01;
end

% Lorenz atractor function (typical values: rho=28; sigma=10; beta=8/3;)
%   x? = -y^2 - z^2 - a*x + a*F 
%   y? = x*y - b*x*z - y + G 
%   z? = b*x*y + x*z -z
f = @(t,X) [...
    -X(2)^2-X(3)^2-a*X(1)+a*F, X(1)*X(2)-b*X(1)*X(3)-X(2)+G, b*X(1)*X(2)+X(1)*X(3)-X(3)];

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