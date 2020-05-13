function hessvec = Adjoint_HessVec( v, Yfull, Lambda )

global N Nt M C nu dt omega

V = reshape(v,N-1,Nt);

%% Solve for w
w=zeros(N-1,Nt);
w(:,1) = zeros(N-1,1); %initial condition
for ii = 2:Nt
    Mat = (1/dt)*M + Ny(Yfull(:,ii)) + nu*C;
    rhs = -(-(1/dt)*M)*w(:,ii-1) - M*V(:,ii);
    
    w(:,ii) = Mat\rhs;
end

%% Solve for p
p=zeros(N-1,Nt);
Mat = (1/dt)*M + Ny(Yfull(:,end)) + nu*C;
rhs = ( dt*M + Nyy(Lambda(:,end)) )*w(:,end);
p(:,end) = Mat'\rhs; %terminal condition

for ii = (Nt-1):-1:1

    Mat = (1/dt)*M + Ny(Yfull(:,ii)) + nu*C;
    rhs = -( -((1/dt)*M) )'*p(:,ii+1) + ( dt*M + Nyy(Lambda(:,ii)) )*w(:,ii);
    
    p(:,ii) = Mat'\rhs;
end

%% Calculate Hess_f*u
HessVec      = zeros(N-1,Nt);
HessVec(:,1) = dt*omega*M*V(:,1); 
for ii = 2:Nt
    HessVec(:,ii) = -M'*p(:,ii) + (dt*omega*M)*V(:,ii);
end

hessvec = reshape(HessVec,(N-1)*Nt,1);

end

