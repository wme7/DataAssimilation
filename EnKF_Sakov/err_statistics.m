function [rmse,m_var,m_std]=err_statistics(nx,nens,nob,H,x,xm,y,yloc)
rmse=0.0;
m_var=0.0;
m_std=0.0;
for p=1:nob
  hx=H(p,:)*x';
  hxm=H(p,:)*xm';
  var=0.0;
  for m=1:nens
    var=var+hx(m)^2;
  end
  var=var/(nens-1);
  rmse=rmse+(y(p)-hxm)^2;
  m_var=m_var+var;
  m_std=m_std+sqrt(var);
end
rmse=sqrt(rmse/nob);
m_var=m_var/nob;
m_std=m_std/nob;
end
