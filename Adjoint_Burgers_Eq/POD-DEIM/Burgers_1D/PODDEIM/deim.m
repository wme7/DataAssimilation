function [P,IND] = deim( U )

[n,m]=size(U);
P=zeros(n,m);
IND=zeros(m,1);

%l=1
[~,ind]=max(abs(U(:,1)));
P(ind,1)=1;
IND(1)=ind;

for l=2:m
   c=(transpose(P(:,1:l-1))*U(:,1:l-1))\(transpose(P(:,1:l-1))*U(:,l)); 
   r=U(:,l)-U(:,1:l-1)*c;
   
   [~,ind]=max(abs(r));

   P(ind,l)=1;
   IND(l)=ind;
end

end