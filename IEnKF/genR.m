function [R] = genR(m)
    
    [R, R1] = qr(randn(m, m));
    for i = 1 : m
        if R1(i, i) < 0
            R(:, i) = -R(:, i);
        end
    end
    
    return
