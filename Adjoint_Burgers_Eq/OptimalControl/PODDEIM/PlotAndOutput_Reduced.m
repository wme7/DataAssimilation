function PlotAndOutput_Reduced(time,space, Yfull, U, Y_red, u, z_red, k, gradf_red, s_red, alpha, count, iter)

global Nt

figure(99)

subplot(2,1,1)
mesh(time, space, [zeros(1,Nt);Yfull;zeros(1,Nt)]);
title('y(t,x)')

subplot(2,1,2)
mesh(time, space, [zeros(1,Nt);U;zeros(1,Nt)]);
title('u(t,x)')

% Some output
if k == 1
    fprintf(1,['     k      f(x_k)    ||gradf(x_k)||   ||s_k||       stepsize    #CG_iter   time(s) \n']);
end
fprintf('  %4d  %12.6e  %12.6e  %12.6e  %7.3f  %8d  %10.4f       - PODDEIM model \n', ...
      k, costFun_Reduced(u,Y_red,z_red), norm(gradf_red), norm(s_red), alpha, iter, count );



end


