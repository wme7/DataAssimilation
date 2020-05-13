function coef=localization(n,x,b)
%covariance localization from Gaspari and Cohn
for i=1:n
  z=min(abs(i-x),n-abs(i-x));
  r=z/b;
  if(z >= 2*b)
      coef(i)=0.0;
  elseif (z >= b && z < 2*b)
      coef(i)=((((r/12.-0.5)*r+0.625)*r+5./3.)*r-5.)*r+4.-2./(3.*r);
  else
      coef(i)=(((-0.25*r+0.5)*r+0.625)*r-5./3.)*r^2+1.;
  end
end
end
      
