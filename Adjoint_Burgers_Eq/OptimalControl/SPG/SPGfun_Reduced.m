function [ f_red,g_red ] = SPGfun_Reduced( u, z_red, y0_red )

global N Nt time space plotOn Upod flag_fg

U = reshape(u,N-1,Nt);


f_red=0;
g_red=zeros((N-1)*Nt,1);

if flag_fg == 1
    Y_red = Burgers_Reduced(U,y0_red);
    [~, g_red] = Adjoint_Grad_Reduced(Y_red, U, z_red);
    f_red = costFun_Reduced(u,Y_red,z_red);
elseif flag_fg == 2
    Y_red = Burgers_Reduced(U,y0_red);
    f_red = costFun_Reduced(u,Y_red,z_red);
elseif flag_fg == 3
    Y_red = Burgers_Reduced(U,y0_red);
    [~, g_red] = Adjoint_Grad_Reduced(Y_red, U, z_red);
end


if plotOn == 1
    figure(111)
    subplot(2,1,1)
    mesh(time, space, [zeros(1,Nt);Upod*Y_red;zeros(1,Nt)]);
    title('y(t,x)')
    
    subplot(2,1,2)
    mesh(time, space, [zeros(1,Nt);U;zeros(1,Nt)]);
    title('u(t,x)')
end

end