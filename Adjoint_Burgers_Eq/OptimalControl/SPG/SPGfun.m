function [ f,g ] = SPGfun( u, z, y0 )

global N Nt time space plotOn flag_fg

U = reshape(u,(N-1),Nt);

f=0;
g=zeros((N-1)*Nt,1);
Lambda = zeros(N-1,Nt); 


if flag_fg == 1
    Yfull = Burgers(U,y0);
    [Lambda, g] = Adjoint_Grad(Yfull, U, z);
    f = costFun(u,Yfull,z);
elseif flag_fg == 2
    Yfull = Burgers(U,y0);
    f = costFun(u,Yfull,z);
elseif flag_fg == 3
    Yfull = Burgers(U,y0);
    [~, g] = Adjoint_Grad(Yfull, U, z);
end    
    
    

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
end

end
