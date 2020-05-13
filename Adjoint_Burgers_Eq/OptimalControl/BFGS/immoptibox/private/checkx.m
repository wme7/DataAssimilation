function  [err, x,n] = checkx(x0)
%CHECKX  Check vector

% Version 04.01.25.  hbn@imm.dtu.dk

err = 0;  sx = size(x0);   n = max(sx);
if  (min(sx) ~= 1) | ~isreal(x0) | any(isnan(x0(:))) | isinf(norm(x0(:))) 
  err = -1;   x = []; 
else
  x = x0(:); 
end