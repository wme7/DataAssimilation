function eval = Ntildey(y_red)

global k1 Bred Upod

tmp = kron(ones(1,k1),Upod*y_red);
eval = Bred*(tmp.*Upod);

end

