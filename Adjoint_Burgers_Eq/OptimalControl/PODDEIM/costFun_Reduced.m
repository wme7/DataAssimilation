function f = costFun_Reduced( u, Y_red, z_red )

global h dt N Nt omega M

U = reshape(u,N-1,Nt);

f = 0;
for t=1:Nt
    f = f + dt*(0.5*Y_red(:,t)'*Y_red(:,t) - h*z_red'*Y_red(:,t) + omega/2*U(:,t)'*M*U(:,t)); 
end

end