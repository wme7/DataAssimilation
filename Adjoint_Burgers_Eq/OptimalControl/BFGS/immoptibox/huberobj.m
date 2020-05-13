function  [f, S, r, J, g] = huberobj(fun, x, gamma, Htype, varargin)
%HUBEROBJ  Value and gradient of the Huber objective function
% for the vector function given by
%            function  [r, J] = fun(x,p1,p2,...)
% p1,p2,... are parameters of the function and J is the Jacobean.
% If nargout < 4, then the function only needs to return r, and the
% gradient is not computed.
%  
% Call    [f, S, r] = huberobj(fun, x, gamma)
%         [f, S, r] = huberobj(fun, x, gamma, Htype, p1,p2,...)
%         [f, S, r, J, g] = huberobj(.....)
%
% Input parameters
% fun   :  Handle to the function.
% x     :  n-vector, argument of fun
% gamma :  Huber threshold.
% Htype :  Choice of Huber function,
%          1 : one-sided,  r_i > 0  are neglected,
%          2 : one-sided,  all r_i <= gamma are active,
%          otherwise,  all abs(r_i) <= gamma  are active (default).
% p1,p2,..  are passed directly to the function FUN .    
%
% Output parameters
% f   :  Huber objective function.
% S   :  Struct with the Huber active set.  Fields
%        s :  Huber sign vector,
%        A :  active set (indices of small elements in r),
%        N :  inactive set (indices of large elements in r),
%        L :  vector  [length(S.A), length(S.N)].
% r,J :  Output from the evaluation of fun .
%        If  nargout < 4 , then  fun  only needs to return  r(x) ,
%        and the gradient is not computed.
% g   :  Gradient of  f .

% Version 10.10.05.  hbn(a)imm.dtu.dk

% Evaluate the function and set rounding error threshold
if  nargout < 4
  r = feval(fun,x, varargin{:});
  tau = 10*eps*norm(r,inf);
else
  [r J] = feval(fun,x, varargin{:});
  tau = 10*eps*max(norm(r,inf), norm(J,inf)*norm(x,inf));
end
if  gamma < tau
  error(sprintf('gamma must be at least %9.2e',tau))
end
thr = gamma + tau;

% Huber sign vector
s = sign(r);
if  nargin < 4 | abs(Htype - 1.5) ~= 0.5,  Htype = 3; end
if  Htype == 1  % Neglect positive contributions
  w = find(-thr <= r & r <= tau);    s(w) = 0;
  p = find(s > 0);   s(p) = 0;
elseif  Htype == 2  % Negative contributions are active
  w = find(r <= thr);   s(w) = 0;
else  % simple Huber function
  w = find(abs(r) <= thr);    s(w) = 0; 
end

% Active set
N = find(s);
S = struct('s',s', 'A',w', 'N',N', 'L',[length(w)  length(N)]);

% Compute function
rA = r(S.A);  rN = r(S.N);  
f = norm(rA)^2/(2*gamma) + norm(rN,1) - .5*gamma*S.L(2);

if  nargout > 3  % gradient with check of rounding errors
  g = J(S.A,:)'*rA/gamma + J(S.N,:)'*s(S.N);
  thr = 10*eps*(norm(J(S.A,:),1)*norm(rA,inf)/gamma + norm(J(S.N,:),1));
  if  norm(g,inf) <= thr,  g = zeros(size(g)); end  
end  % gradient