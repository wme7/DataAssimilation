function PlotAndOutput(time,space, Yfull, z, Lambda, U, k, gradf, s, alpha, count, iter)

global Nt N plotOn

if plotOn == 1
    if k==0
        figure(66)

        subplot(2,2,1)
        mesh(time, space, [zeros(1,Nt);Yfull;zeros(1,Nt)]);
        title('y(t,x)')

        subplot(2,2,2)
        mesh(time, space, kron(ones(1,Nt),[0;z;0]));
        title('z(t,x)')

        subplot(2,2,4)
        mesh(time, space, [zeros(1,Nt);U;zeros(1,Nt)]);
        title('u(t,x)')
    else
        figure(66)
        subplot(2,2,1)
        mesh(time, space, [zeros(1,Nt);Yfull;zeros(1,Nt)]);
        %zlim([0 1.5])
        title('y(t,x)')

        subplot(2,2,2)
        mesh(time, space, kron(ones(1,Nt),[0;z;0]));
        %zlim([0 1.5])
        title('z(t,x)')

        subplot(2,2,3)
        mesh(time, space, [zeros(1,Nt);Lambda;zeros(1,Nt)]);
        title('\lambda(t,x)')
        subplot(2,2,4)

        mesh(time, space, [zeros(1,Nt);U;zeros(1,Nt)]);
        title('u(t,x)')

        % Some output
        if k == 1
            fprintf(1,['     k      f(x_k)    ||gradf(x_k)||   ||s_k||     stepsize    #CG_iter   time(s)\n']);
        end
        fprintf('  %4d  %12.6e  %12.6e  %12.6e  %7.3e  %8d %11.4f     - full model \n', ...
              k, costFun(reshape(U,(N-1)*Nt,1),Yfull,z), norm(gradf), norm(s), alpha, iter, count );

    end
end

end




