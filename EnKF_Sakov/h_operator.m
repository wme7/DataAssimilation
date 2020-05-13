function H=h_operator(N,P,obloc)
H=zeros(P,N);
for p=1:P
for i=1:N-1
 if(obloc(p)>=i && obloc(p)<=i+1)
   H(p,i)=1-obloc(p)+i;
   H(p,i+1)=obloc(p)-i;
 end
end

end
