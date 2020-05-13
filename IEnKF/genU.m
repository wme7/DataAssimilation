function [U] = genU(m)
  
    v = ones(m, 1);
    [V, TMP] = svd(v);
    Rm1 = genR(m - 1);
    U = V * blkdiag(1, Rm1) * V';
  
    return
