function J = JacFN(t,y)

global A N

fprime = @(v) ((v-0.1).*(1-v) + v.*(1-v) - v.*(v-0.1));

Fprime = [fprime(y(1:N));zeros(N,1)];

J = A + spdiags(Fprime,0,2*N,2*N);

end