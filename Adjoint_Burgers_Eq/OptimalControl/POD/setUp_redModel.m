function [z_red, y0_red] = setUp_redModel(Yfull, z, y0)

global M B C k1 plotOn N Nt Mred1 Bred Cred Upod TMP

% SVD of snapshot matrices + DEIM
R = chol(M);
[Uy,Sy,~]=svd(R*Yfull);
Upod=R\Uy(:,1:k1);


if plotOn ==1
   figure(55)
   index=1:min(N-1,Nt);
   semilogy(index,diag(Sy(index,index)),'x'); %investigate basis dimensions
end


% Pre-compute matrices
Mred1 = Upod'*M;
Cred  = Upod'*C*Upod;
Bred  = Upod'*B;
TMP = zeros(k1,k1,k1);
for ii=1:k1
    for jj=1:ii
        TMP(:,ii,jj)=Bred*(Upod(:,ii).*Upod(:,jj));      
    end
end
z_red = Upod'*z; %just because z is multiplied by y in J(.)
         
% obtain reduced variables
y0_red = Upod'*M*y0;

end

