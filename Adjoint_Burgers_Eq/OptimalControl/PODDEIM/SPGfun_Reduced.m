function [ f_red,g_red ] = SPGfun_Reduced( u_red, z_red, y0_red )

global k3 Nt time space plotOn Upod Uupod

U_red = reshape(u_red,k3,Nt);

Y_red = Burgers_Reduced(U_red,y0_red);
[~, g_red] = Adjoint_Grad_Reduced(Y_red, U_red, z_red);
f_red = costFun_Reduced(u_red,Y_red,z_red);

if plotOn == 1
    figure(111)
    subplot(2,1,1)
    mesh(time, space, [zeros(1,Nt);Upod*Y_red;zeros(1,Nt)]);
    title('y(t,x)')
    
    subplot(2,1,2)
    mesh(time, space, [zeros(1,Nt);Uupod*U_red;zeros(1,Nt)]);
    title('u(t,x)')
end

end