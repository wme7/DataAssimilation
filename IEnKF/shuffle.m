function [pack] = shuffle(n)
    pack = 1 : n;
    for i = 1 : n
        ii = floor(rand * n) + 1;
        val_ii = pack(ii);
        pack(ii) = pack(i);
        pack(i) = val_ii;
    end
  
    return
