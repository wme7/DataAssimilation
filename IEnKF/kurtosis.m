function[k] = kurtosis(X)
    
    if length(X) == 0
        k = NaN;
    elseif sum(size(X) > 1) == 1
        a = X - mean(X);
        k = sum(a .^ 4) / sum(a .^ 2) ^ 2 * length(X);
    else
        A = X - repmat(mean(X), size(X, 1), 1);
        k = sum(A .^ 4) ./ sum(A .^ 2) .^ 2 * size(X, 1);
    end
    
    return
