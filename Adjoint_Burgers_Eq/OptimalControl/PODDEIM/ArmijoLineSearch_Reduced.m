function [alpha, Y, u] = ArmijoLineSearch_Reduced( f, Y, y0,  z, u0, s, sTgradf )          
%
%  Inputs:
%  f  - function handle that should be minimized
%  u0 - starting point
%  s - search direction
%  sTgradf - precomputed < s , grad_f(u0) >
%
%  Output: alpha, the optimal step size in direction s
%


global tol_line min_alpha N Nt


%% Init
alpha = 1;
f0    = f(u0, Y, z);

u = u0+alpha*s;
Y = Burgers_Reduced(reshape(u,N-1,Nt),y0);

feval = f(u, Y, z);


%% Start linear search loop
while (feval > f0 + alpha*tol_line*sTgradf) & alpha > min_alpha 
    
    % reduce step size 
    alpha  = alpha/2;
    
    % Calc new u, and corr. Y 
    u = u0+alpha*s;
    Y = Burgers_Reduced(reshape(u,N-1,Nt),y0);
    feval = f(u, Y, z);

   
end


end
