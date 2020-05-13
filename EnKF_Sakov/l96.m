function x=l96(x,dt,N,a,b,d,F)
smalldt=0.05;
for t=1:dt/smalldt
  for i=1:N
    x1(:,:,i)=x(:,:,i)+(x(:,:,rema(i-1,N)).*(a*x(:,:,rema(i+1,N))-b*x(:,:,rema(i-2,N)))-d*x(:,:,i)+F)*smalldt/2;
  end
  for i=1:N
    x2(:,:,i)=x(:,:,i)+(x1(:,:,rema(i-1,N)).*(a*x1(:,:,rema(i+1,N))-b*x1(:,:,rema(i-2,N)))-d*x1(:,:,i)+F)*smalldt/2;
  end
  for i=1:N
    x3(:,:,i)=x(:,:,i)+(x2(:,:,rema(i-1,N)).*(a*x2(:,:,rema(i+1,N))-b*x2(:,:,rema(i-2,N)))-d*x2(:,:,i)+F)*smalldt;
  end
  for i=1:N
    x4(:,:,i)=x(:,:,i)+(x(:,:,rema(i-1,N)).*(a*x(:,:,rema(i+1,N))-b*x(:,:,rema(i-2,N)))-d*x(:,:,i)+F)*smalldt/6+...
                       (x1(:,:,rema(i-1,N)).*(a*x1(:,:,rema(i+1,N))-b*x1(:,:,rema(i-2,N)))-d*x1(:,:,i)+F)*smalldt/3+...
                       (x2(:,:,rema(i-1,N)).*(a*x2(:,:,rema(i+1,N))-b*x2(:,:,rema(i-2,N)))-d*x2(:,:,i)+F)*smalldt/3+...
                       (x3(:,:,rema(i-1,N)).*(a*x3(:,:,rema(i+1,N))-b*x3(:,:,rema(i-2,N)))-d*x3(:,:,i)+F)*smalldt/6;
  end
  x=x4;
end

end

function r=rema(i,N)
  r=mod(i-1,N)+1;
end
