function [Q,R,status]=inflation_additive(N,M,P,H,uf,uf1,ua,y,y1,Q,R)

    tau=10;

    xfm=mean(squeeze(uf),1);
    xf1m=mean(squeeze(uf1),1);
    xam=mean(squeeze(ua),1);
    xf=squeeze(uf)-repmat(xfm,M,1);
    xf1=squeeze(uf1)-repmat(xf1m,M,1);
    xa=squeeze(ua)-repmat(xam,M,1);

    e=y-xfm*H';
    e1=y1-xf1m*H';
    Ke=xam-xfm;

    if(cond(xa)>1e10)
      status=1;
      return
    end

    HF=H*xf1'/xa';
    if(cond(HF)>1e10)
      status=1;
      return
    end

    Q1=(HF\e1'*e+Ke'*e)/H'-xf1'*xf1/N;
    Q=Q+(Q1-Q)/tau;
    Q=(Q+Q')/2;
    Q=real(Q);

%    R1=e'*e-H*xf1'*xf1*H'/N;
%    R=R+(R1-R)/tau;
%    R=real(R);

    status=0;

end
