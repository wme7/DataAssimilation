function [z_red, y0_red] = setUp_redModel(Yfull, z, y0)

global M B C k1 k2 plotOn N Nt Mred1 Bred Cred Fred Upod TMP

% SVD of snapshot matrices + DEIM
R = chol(M);
[Uy,Sy,~]=svd(R*Yfull);
Upod=R\Uy(:,1:k1);
% Upod = setupFEMPOD(R*Yfull);

[Uf,Sf,~]=svd(Yfull.^2);
Udeim=Uf(:,1:k2);
[~,ind]=deim(Udeim);

if plotOn ==1
   figure(55)
   index=1:min(N-1,Nt);
   semilogy(index,diag(Sy(index,index)),'x',index,diag(Sf(index,index)),'rx'); %investigate basis dimensions
end

% Pre-compute matrices
Mred1 = Upod'*M;
% Bred  = Upod'*B*Udeim*inv(Pdeim'*Udeim);
Bred  = Upod'*B*Udeim/Udeim(ind,:);
Cred  = Upod'*C*Upod;
% Fred  = Pdeim'*Upod;
Fred = Upod(ind,:);
TMP = zeros(k1,k1,k1);
for ii=1:k1
    for jj=1:ii
        TMP(:,ii,jj)=Bred*(Fred(:,ii).*Fred(:,jj));      
    end
end
z_red = Upod'*z; %just because z is multiplied by y in J(.)

% obtain reduced variables
y0_red = Upod'*M*y0;
% Y_red = Burgers_Reduced(U_red,y0_red);

end

