function [ f,g ] = BFGSfun( u, z, y0 )

global N Nt iter_tmp time space plotOn

U = reshape(u,(N-1),Nt);

Yfull = Burgers(U,y0);
[Lambda, g] = Adjoint_Grad(Yfull, U, z);
f = costFun(u,Yfull,z);

if plotOn == 1
    figure(1)
    subplot(2,2,1)
    mesh(time, space, [zeros(1,Nt);Yfull;zeros(1,Nt)]);
    zlim([0 1.5])
    title('y(t,x)')
    subplot(2,2,2)
    mesh(time, space, kron(ones(1,Nt),[0;z;0]));
    zlim([0 1.5])
    title('z(t,x)')
    subplot(2,2,3)
    mesh(time, space, [zeros(1,Nt);Lambda;zeros(1,Nt)]);
    title('\lambda(t,x)')
    subplot(2,2,4)
    mesh(time, space, [zeros(1,Nt);U;zeros(1,Nt)]);
    title('u(t,x)')

    if iter_tmp == 0
        fprintf(1,['     #     f(x_k)    ||gradf(x_k)||  \n']);
    end
    
    fprintf('  %4d  %12.6e  %12.6e \n', ...
          iter_tmp, f , norm(g));
    
end
 
if iter_tmp == 1
    u_ss=u;
    save ./BFGS/snapshots.mat u_ss
end


iter_tmp = iter_tmp +1;

end
