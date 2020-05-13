function eval = rhsPODDEIM(t,a)

global A C N eps dx E U PP

y=U*a;

ifun = @(t) (50000*t.^3 .* exp(-15*t));
ffun = @(v) (v.*(v-0.1).*(1-v)); %the nonlin function

f=ffun(y(1:N));
Fdeim=PP*f;

F = [Fdeim;zeros(N,1)]; %non-linearity

g = zeros(2*N,1);
g(1,1) = eps^2/dx * ifun(t); %set Neumann BC

eval = U'*inv(E)*(A*y + C + g + F);

end