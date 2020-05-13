function eval = Ntildey(y_red)

global k1 Bred Fred

tmp = kron(ones(1,k1),Fred*y_red);
eval = Bred*(tmp.*Fred);

end

