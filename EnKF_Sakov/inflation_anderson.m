function [inf,var_inf]=inflation_anderson(nx,nens,nob,H,inf,var_inf,obserr,x,xm,y,yloc,localize,l_coef)

pi=4.0*atan(1.0);

for p=1:nob
  hx=H(p,:)*x';
  hxm=H(p,:)*xm';

  dist=abs(hxm-y(p));

  var=0.0;
  cov=0.0;
  for m=1:nens
    var=var+hx(m)*hx(m);
    cov=cov+x(m,:)*hx(m);
  end
  var=var/(nens-1);
  cov=cov/(nens-1);

  if(localize)
    covloc=localization(nx,yloc(p),l_coef);
  else
    covloc=1.0;
  end

  gama=covloc.*cov/var;
  inf_o=(1+gama.*(sqrt(inf)-1)).^2;
  th=sqrt(inf_o*var+obserr^2);
  lm=exp(-0.5*(dist./th).^2)./(sqrt(2.*pi)*th);
  lp=lm.*((dist./th).^2-1)./th;
  lp=lp*0.5*var.*gama.*(1-gama+gama.*sqrt(inf))./(th.*sqrt(inf));
  bb=lm./lp-2*inf;
  cc=inf.^2-var_inf-lm.*inf./lp;
  inf1=0.5*(-bb+sqrt(bb.^2-4*cc));
  inf2=0.5*(-bb-sqrt(bb.^2-4*cc));
  for i=1:nx
    if(gama(i)>0.0)
      if(abs(inf1(i)-inf(i)) > abs(inf2(i)-inf(i)))
        if(inf2(i)>0)
            inf(i)=inf2(i);
        end
      else
        if(inf1(i)>0)
            inf(i)=inf1(i);
        end
      end
    end
  end

end

%var_b=diag(x'*x/(nens-1));
%sig1=sqrt((inf+sqrt(var_inf)).*var_b'+obserr^2);
%sig2=sqrt(inf.*var_b'+obserr^2);
%n1=(1./sig1).*exp(-0.5./(sig1.^2)-0.5);
%n2=(1./sig2).*exp(-0.5./(sig2.^2));
%ratio=n1./n2;
%var_inf=-var_inf./(2.*log(ratio));
%var_inf=1.1*var_inf; %inflate var_inf to avoid filter divergence in inf;

end
