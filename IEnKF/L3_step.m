function [x] = L3_step(x)
  
    x = rk4step(@L3, 0, 0.01, x);
  
    return
