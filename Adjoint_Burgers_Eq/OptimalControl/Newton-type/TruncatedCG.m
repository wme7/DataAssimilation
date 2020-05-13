function [sk,i] = TruncatedCG( Hess, s0, rhs, max_it ) 

global eta

sk = s0;
p0 = rhs;
r0 = rhs;

rold = r0;
p = p0;

for i=1:max_it

    if norm(rold) < eta*norm(r0)
        if i == 1
            sk = rhs;
        end
        
        return;
    end
    
    q = Hess(p);
    
    if transpose(p)*q < 0
        if i == 1
            sk = rhs;
        end
        
        return;

    end
    
    gamma = norm(rold)^2/(transpose(p)*q);
    
    sk = sk + gamma*p; 
    rnew = rold - gamma*q;
    beta = norm(rnew)^2/norm(rold)^2;
    p = rnew + beta*p;
    
    rold = rnew;
end




end