function X = Lorenz96(N,T,dt,F,EPS)
% Lorenz96 Function: solves explicitly the lorenz 96 system of PDE's
%
%       X = Lorenz96(T,dt,F,epsilon)
%
%       where: 
%
%       N       - Number of variables
%       T       - time range [0,tEnd]
%       dt      - time step
%       F       - Forcing term
%       EPS     - perturbation
%       X       - output vectors
%
% Example.
%        [x y z] = Lorenz96(36, 1.0, 0.01, 8.0, 0.01); or
%        [x y z] = Lorenz96();
%        plot3(x y z);
%
if nargin == 0
    N    = 36;
    T    = [0,30];
    dt   = 0.01;
    F    = 8.0;
    EPS  = 0.01;
end

X0 = F * ones(N,1);     % Initial state (equilibrium)
X0(20) = X0(20) + EPS;  % Add small perturbation to 20th variable
t = (T(1):dt:T(end))';	% build time domain
dimX = length(X0);      % spatial size
dimT = length(t);       % temporal size
X = zeros(dimT,dimX);   % solution array 
X(1,:) = X0;            % Load initial condition

% Produce an explicit RK4 time integration
for n = 1:dimT-1
    k1 = f(t(n,:)       ,X(n,:)          ,F,N);
    k2 = f(t(n,:)+0.5*dt,X(n,:)+0.5*dt*k1,F,N);
    k3 = f(t(n,:)+0.5*dt,X(n,:)+0.5*dt*k2,F,N);
    k4 = f(t(n,:)+  dt  ,X(n,:)+  dt*k3  ,F,N);
    X(n+1,:) = X(n,:) + dt*(k1+2*k2+2*k3+k4)/6;
end

end

function d = f(~,X,F,N)
    % organize as row vectors
    d = zeros(1,N);
    % First, compute the 3 problematic cases: 1,2,N
    d(1)=(X(2)-X(N-1))*X(N)-X(1);
    d(2)=(X(3)-X(N))*X(1)-X(2);
    d(N)=(X(1)-X(N-2))*X(N-1)-X(N);
    % then the general cases
    for i=3:N-1, d(i)=(X(i+1)-X(i-2)).*X(i-1)-X(i); end
    %Lastly add the Forcing term
    d = d + F;
end