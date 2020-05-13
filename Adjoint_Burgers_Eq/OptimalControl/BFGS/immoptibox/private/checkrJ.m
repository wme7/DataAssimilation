function  [err, f, r, J] = checkrJ(fun, x, varargin)
%CHECKRJ  Check Matlab function which is called by a 
% nonlinear least squares solver.

% Version 10.09.20.  hbn(a)imm.dtu.dk

err = 0;   f = NaN;  n = length(x);
if  nargout > 3    % Check  r  and  J
  [r J] = feval(fun,x,varargin{:});
  sr = size(r);   sJ = size(J);
  if  sr(2) ~= 1 | ~isreal(r) | any(isnan(r(:))) | any(isinf(r(:)))
    err = -2;  return, end
  if  ~isreal(J) | any(isnan(J(:))) | any(isinf(J(:)))
    err = -3;  return, end
  if  sJ(1) ~= sr(1) | sJ(2) ~= n
    err = -4;  return, end
  
else  % only check  r
  r = feval(fun,x,varargin{:});
  sr = size(r);   
  if  sr(2) ~= 1 | ~isreal(r) | any(isnan(r(:))) | any(isinf(r(:)))
    err = -2;  return, end
end

% Objective function
f = (r' * r)/2;
if  isinf(f),  err = -5; end