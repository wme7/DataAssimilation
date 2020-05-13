function [inf,var_inf]=inflation_anderson_old(inf,var_inf,obserr,var_b,dist)
a=1.0;
b=-obserr*obserr-inf*var_b;
c=0.5*var_inf*var_b*var_b;
d=-0.5*var_inf*var_b*var_b*dist*dist;
p=-b/(3.0*a);
q=p^3+(b*c-3.0*a*d)/(6*a*a);
r=c/(3.0*a);
x=(q+(q*q+(r-p*p)^3)^(1./2.))^(1./3.)+(q-(q*q+(r-p*p)^3)^(1./2.))^(1./3.)+p;
inf=(x-obserr^2)/var_b;
%constant var_inf, or use the following to update var_inf
%sig1=sqrt((inf+sqrt(var_inf))*var_b+obserr^2);
%sig2=sqrt(inf*var_b+obserr^2);
%n1=(1./sig1)*exp(-0.5/(sig1^2)-0.5);
%n2=(1./sig2)*exp(-0.5/(sig2^2));
%ratio=n1/n2;
%var_inf=-var_inf/(2.*log(ratio));
%var_inf=1.2*var_inf %inflate var_inf to avoid filter divergence in inf;
end
