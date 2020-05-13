function [y] = rmse(x)
    [n, m] = size(x);
    
    if n == 1 | m == 1
        y = sqrt(sum(x .^ 2) / length(x));
    else
        y = zeros(m, 1);
        for i = 1 : m
            y(i) = sqrt(sum(x(:, i) .^ 2) / n);
        end
    end
    
    return
