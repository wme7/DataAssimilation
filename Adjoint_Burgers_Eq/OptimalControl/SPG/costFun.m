function f = costFun( u, Y, z )

global dt Nt omega M h N

U = reshape(u,N-1,Nt);

f = 0;
for t=1:Nt
    f = f + dt*(0.5*Y(:,t)'*M*Y(:,t) - h*z'*Y(:,t) + omega/2*U(:,t)'*M*U(:,t)); 
end


