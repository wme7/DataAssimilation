function Yfull = Burgers(y0)

global tol_newton Nt N max_newton

rhs = zeros(N-1,1); % zero rhs in Burgers eqn

Yfull      = zeros(N-1,Nt);
Yfull(:,1) = y0;

for tt=2:Nt %time loop
    y_tmp1=Yfull(:,tt-1); %nice guess  
    err=1; 
    iter = 1;
    while (err > tol_newton && iter < max_newton)%inner Newton method
       f=fNewton(Yfull(:,tt-1),y_tmp1,rhs); 
       J=JacNewton(y_tmp1);

       ytilde=J\-f;
       y_tmp2=y_tmp1+ytilde;

       err=norm(y_tmp2-y_tmp1)/norm(y_tmp2);
       y_tmp1=y_tmp2;  
       
       iter = iter+1;
    end
    Yfull(:,tt)=y_tmp2;
end

end
