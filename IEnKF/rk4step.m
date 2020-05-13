function [yend] = rk4step(f, x, xend, y)
  
    h = xend - x;
  
    k1 = h * feval(f, x, y);
    k2 = h * feval(f, x + 0.5 * h, y + 0.5 * k1);
    k3 = h * feval(f, x + 0.5 * h, y + 0.5 * k2);
    k4 = h * feval(f, x + h, y + k3);

    yend = y + (k1 + 2 * (k2 + k3) + k4) / 6;
  
    return
