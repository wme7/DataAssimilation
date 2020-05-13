function [posterior,lambda,r_coef,cr,infx]=enkf(prior,N,P,M,H,ob,obloc,R,Q,inf_additive,inf_type,inf,inf_var,relax_type,lambda0,r_coef0,localize,l_coef,tau)

x=squeeze(prior);
y=squeeze(ob);
xm=mean(x,1);
x=x-repmat(xm,M,1);

x_orig=x;

std_b=repmat(std(x),M,1);
omb=y-xm*H';
HPfH=H*x'*x*H'/(M-1);

%%%%%%%%%%%%Inflation%%%%%%%%%%%
infx=1.0;
if(strcmp(inf_type,'regular'))
  x=sqrt(inf)*x;
  infx=sqrt(inf);
end
if(strcmp(inf_type,'anderson')) 
  [inf,inf_var]=inflation_anderson(N,M,P,H,inf,inf_var,mean(sqrt(diag(R))),x,xm,y,obloc,localize,l_coef);
  for m=1:M
    x(m,:)=sqrt(inf).*x(m,:);
  end
  infx=sqrt(inf);
end

if(inf_additive)
  [u,l] = eig(Q);
  sqrtinf = real(u*diag(sqrt(max(diag(l),0)))*inv(u));
  x=x+randn(M,N)*sqrtinf';
end

if(strcmp(inf_type,'new1'))
  x=x*sqrt(al);
end

%%%%%%%%%%EnKF%%%%%%%%%%%%%%%%
%make sure obserr has zero mean
%y=y-mean(y);

for p=1:P
  obserr=sqrt(squeeze(R(p,p)));
  hx=H(p,:)*x';
  hxm=H(p,:)*xm';
  var=hx*hx'/(M-1);
  d=var+obserr^2;
  alpha=1.0/(1.0+sqrt(obserr^2/d));
  if(localize)
    rho=localization(N,obloc(p),l_coef);
  else
    rho=ones(1,N);
  end
  K=x'*hx'/(M-1)/d;
  K=K.*rho';
  x=x-alpha*hx'*K';
  xm=xm+(y(p)-hxm)*K';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

std_a=repmat(std(x),M,1);
oma=y-xm*H';
amb=omb-oma;
HPaH=H*x'*x*H'/(M-1);
cr=sqrt((trace(HPfH)+trace(R))/(omb*omb'));
%cr=sqrt(trace(HPaH)/((omb-oma)*oma'));

%%%%%%%%%%Relaxation%%%%%%%%%%
r_coef=r_coef0;
lambda=lambda0;
if(strcmp(relax_type,'rtps'))
  x=x.*(r_coef.*(std_b-std_a)./std_a+1);
  infx=(r_coef.*(std_b(1,:)-std_a(1,:))./std_a(1,:)+1);
end
if(strcmp(relax_type,'rtpp'))
  x=(1-r_coef)*x+r_coef*x_orig;
end
if(strcmp(relax_type,'rtps_adaptive'))
%  sig_ratio=sqrt(trace((HPfH-HPaH)./HPaH)/P);
  sig_ratio=(sqrt(trace(HPfH)/P)-sqrt(trace(HPaH)/P))/sqrt(trace(HPaH)/P);
%  lambda=sqrt((omb*omb'-trace(R))/trace(HPfH));  %match prior spread
  lambda=sqrt((omb-oma)*oma'/trace(HPaH));       %match posterior spread
  lambda=lambda0+(lambda-lambda0)/tau;
  r_coef=(lambda-1)/sig_ratio;
%  r_coef=min(r_coef,1);
%  r_coef=r_coef0+(r_coef-r_coef0)/tau;
  x=x.*(r_coef.*(std_b-std_a)./std_a+1);
  infx=(r_coef.*(std_b(1,:)-std_a(1,:))./std_a(1,:)+1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

posterior(1,:,:)=x+repmat(xm,M,1);

end
