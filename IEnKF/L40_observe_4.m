function [y] = L40_observe(x, r)

    xx = x(1 : 4 : 40);
    if r == 0
        y = xx;
    else
        y = xx + r * randn(size(xx));
    end
    
    return
