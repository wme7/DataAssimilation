function  [err, f, g, H] = checkfgH(fun, x, varargin)
%CHECKFGH  Check Matlab function which is called by a 
% general optimization function

% Version 09.10.02.  hbn(a)imm.dtu.dk

err = 0;  
if  nargout > 3
  [f g H] = feval(fun, x, varargin{:});
elseif  nargout > 2
  [f g] = feval(fun, x, varargin{:});  H = [];
else
  f = feval(fun, x, varargin{:});  g = [];  H = [];
end

sf = size(f);     
if  any(sf ~= 1) | ~isreal(f) | any(isnan(f(:))) | any(isinf(f(:)))
  err = -2;  return, end  % f is not a real valued scalar

if  ~isempty(g)
  sg = size(g);
  if  ~isreal(g) | any(isnan(g(:))) | any(isinf(g(:))) | min(sg) ~= 1
    err = -3;  return, end  % g is not a real valued vector
  if  max(sg) ~= length(x)
    err = -4;  return, end  % dimension mismatch in x, g
end

if  ~isempty(H)
  sH = size(H);
  if  ~isreal(H) | any(isnan(H(:))) | any(isinf(H(:)))
    err = -3;  return, end  % H is not a real valued matrix
  if  any(sH ~= length(x))
    err = -4;  return, end  % dimension mismatch in x, H
  E = H - H';
  if  norm(E(:)) > 10*eps*norm(H(:))
    err = -5;  return,  end  % H is not symmetric
  H = (H + H')/2;
end