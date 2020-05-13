function [y] = L40_observe(x, r)

    if r == 0
        y = x;
    else
        y = x + r * randn(size(x));
    end
    
    return
