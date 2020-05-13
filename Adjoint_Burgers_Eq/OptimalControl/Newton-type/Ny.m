function eval = Ny( y )

global B N

eval = B.*kron(ones(1,N-1),y)';


end

