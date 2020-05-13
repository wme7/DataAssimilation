function [Lambda, gradf] = Adjoint_Grad(Yfull, U, z)


global N Nt dt M B C nu h omega


%% Compute the adjoint Lambda
Lambda = zeros(N-1,Nt); 
% terminal condition
Mat = (1/dt)*M + Ny(Yfull(:,end)) + nu*C ;
rhs = -dt*( M*Yfull(:,end) - h*z );
Lambda(:,end) = Mat'\rhs;

% Solving backwards
for ii=(Nt-1):-1:1
    Mat = (1/dt)*M + Ny(Yfull(:,ii)) + nu*C;
    rhs = -( (-(1/dt)*M)'*Lambda(:,ii+1) ) - dt*( M*Yfull(:,ii) - h*z );
    
    Lambda(:,ii) = Mat'\rhs;
end


%% Compute gradf
Gradf = zeros(N-1,Nt); 
Gradf(:,1) = dt*omega*M*U(:,1);
for ii=2:Nt
    Gradf(:,ii) = dt*omega*M*U(:,ii) - M'*Lambda(:,ii);
end

gradf  = reshape(Gradf,(N-1)*Nt,1);

end