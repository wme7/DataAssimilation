function eval = rhs(t,y)

global A C N eps dx

ifun = @(t) (50000*t.^3 .* exp(-15*t));
ffun = @(v) (v.*(v-0.1).*(1-v));

F = [ffun(y(1:N));zeros(N,1)]; %non-linearity

g = zeros(2*N,1);
g(1,1) = eps^2/dx * ifun(t); %set Neumann BC

eval = A*y + C + g + F;

end