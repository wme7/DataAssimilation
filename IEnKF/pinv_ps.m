function [B] = pinv(A)
    
    DEBUG = 0;
    
    [U, L, V] = svd(A);
    
    if size(A, 2) == 1
        l = L;
    else
        l = diag(L);
    end
    valid = find(l / l(1) > 1e-10);
    invalid = find(l / l(1) <= 1e-10);
    
    m = size(A, 2);
    for i = valid'
        L(i, i) = 1 / l(i);
    end
    for i = invalid'
        L(i, i) = 0;
    end
    
    if DEBUG
        fprintf('  U = %d x %d, L = %d x %d, V = %d x %d\n', size(U, 1), size(U, 2), size(L, 1), size(L, 2), size(V, 1), size(V, 2));
    end

    B = V * L' * U';
    
    return
    