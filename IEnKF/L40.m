function [y1] = L40(t, y)

    n = length(y);
    y1 = zeros(n, 1);
    yp1 = zeros(n, 1);
    ym1 = zeros(n, 1);
    ym2 = zeros(n, 1);
  
    yp1(1 : n - 1) = y(2 : n);
    yp1(n) = y(1);
    ym1(2 : n) = y(1 : n - 1);
    ym1(1) = y(n);
    ym2(3 : n) = y(1 : n - 2);
    ym2(1 : 2) = y(n - 1 : n);
  
    y1 = (yp1 - ym2) .* ym1 - y + 8.0;
  
    return
