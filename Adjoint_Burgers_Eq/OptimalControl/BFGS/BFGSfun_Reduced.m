function [ f_red,g_red ] = BFGSfun_Reduced( u, z_red, y0_red )

global N Nt iter_tmp plotOn time space Upod

U = reshape(u,N-1,Nt);

Y_red = Burgers_Reduced(U,y0_red);
[~, g_red] = Adjoint_Grad_Reduced(Y_red, U, z_red);
f_red = costFun_Reduced(u,Y_red,z_red);


if plotOn == 1
    figure(2)
    subplot(2,1,1)
    mesh(time, space, [zeros(1,Nt);Upod*Y_red;zeros(1,Nt)]);
    title('y(t,x)')
    
    subplot(2,1,2)
    mesh(time, space, [zeros(1,Nt);U;zeros(1,Nt)]);
    title('u(t,x)')

    if iter_tmp == 0
        fprintf(1,['     #     f(x_k)    ||gradf(x_k)||  \n']);
    end
    
    fprintf('  %4d  %12.6e  %12.6e \n', ...
          iter_tmp, f_red , norm(g_red));
    
end
iter_tmp = iter_tmp +1;

end

