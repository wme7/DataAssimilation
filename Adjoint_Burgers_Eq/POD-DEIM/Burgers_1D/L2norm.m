function eval = L2norm( F, G )

global h dt Nt N

eval = 0;
for t=1:Nt
    for i=1:N-1
        eval = eval + dt*h*( F(i,t) - G(i,t) )^2;
    end
end
eval = sqrt(eval);
end

