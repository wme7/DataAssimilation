function  [X, S, info, perf] = nonlinhuber(fun, x0, opts, varargin)
%NONLINHUER  Levenberg-Marquardt type method for Huber estimation.
% Find  xm = argmin{f_gamma(x)} , where  x  is an n-vector and
% f_gamma(x) = 0.5/gamma * sum(phi(r_i(x))) , the Huber objective function.
% The functions  r_i(x) (i=1,...,m) and the Jacobian matrix  J(x)  
% (with elements  J(i,j) = Df_i/Dx_j ) must be given by a MATLAB function 
% with declaration
%            function  [r, J] = fun(x,p1,p2,...)
% p1,p2,... are parameters of the function.  
%
% Call  [X, S, info] = nonlinhuber(fun, x0)
%       [X, S, info] = nonlinhuber(fun, x0, opts, p1,p2,...)
%       [X, S, info, perf] = nonlinhuber(.....)
%
% Input parameters
% fun  :  Handle to the function.
% x0   :  Starting guess for Huber estimator.
% opts :  Either a struct with fields 'gamma', 'Htype', 'tau', 'tolg', 
%         'tolx' and 'maxeval', or a vector with the values of these 
%         options, opts = [gamma  Htype  tau  tolg  tolx  maxeval].
%         gamma:  Huber threshold.
%         Htype:  Choice of Huber function,
%                 Htype = 1: one-sided,  r_i > 0  are neglected,
%                         2: one-sided,  all r_i <= gamma are active,
%                         otherwise,  all abs(r_i) <= gamma  are active.
%           tau:  used to get starting value for Marquardt parameter: 
%                     mu = tau * max( [J(x0)'*J(x0)]_(ii) ) .
%         The remaining options are used in stopping criteria:
%             ||g||_inf <= tolg                            or 
%             ||dx||_2  <= tolx*(tolx + ||x||_2)           or
%             no. of function evaluations exceeds  maxeval .
%         Default  gamma = ||r(x0)||inf,  Htype = 0,  
%         tau = 1e-3,  tolg = 1e-4,  tolx = 1e-8, maxeval = 100.  
%         If the input opts has less than 6 elements, it is augmented by 
%         the default values.  Also, zeros and negative elements are 
%         replaced by the default values.
% p1,p2,..  are passed directly to the function FUN .    
%
% Output parameters
% X   :  If perf is present, then array, holding the iterates columnwise, 
%        otherwise computed solution. 
% S   :  Struct with the Huber active set for the last col. in X.  Fields
%        s:  Huber sign vector,
%        A:  active set (indices of small components in r),
%        N:  indices of large components,
%        L:  vector [length(S.A), length(S.N)].
% info:  Performance information, vector with 7 elements:
%        info(1:4) = final values of 
%             [f(x)  ||g(x)||inf  ||dx||2  mu/max( [J(x)'*J(x)]_(ii) )] ,
%        info(5:6) = no. of iteration steps and function evaluations.
%        info(7) = 1 :  Stopped by small gradient
%                  2 :  Stopped by small x-step
%                  3 :  No. of function evaluations exceeded.
% perf :  Struct with fields
%          f :  values of  f(xk) ,
%         ng :  values of  || g(xk) ||_inf ,
%         mu :  values of Marquardt damping parameter.

% Version 10.11.16.  hbn(a)imm.dtu.dk

% Check options
if  nargin < 3,  opts = []; end
opts  = checkopts('nonlinhuber', opts);  % use default options where required
gamma = opts(1);  Htype = opts(2);
tau = opts(3);    tolg = opts(4);  tolx = opts(5);  maxeval = opts(6);

%  Initialize 
x = x0(:);  
if  gamma == 0  % get gamma such that all eqs are active
  gamma = norm(feval(fun,x, varargin{:}), inf);
  neval = 1;
else,  neval = 0; end 
[f S r J g] = huberobj(fun,x,gamma,Htype,varargin{:});  neval = neval + 1;
A = J'*J;  mu = tau * max(diag(A));
ng = norm(g,inf);
nit = 1;  % 1 + no. of iterations   
nu = 2;   nh = 0;  stop = 0;  
Trace = nargout > 3;
if  Trace
  o = ones(1,maxeval);
  X = x*o;  perf = [f; ng; mu]*o;
end 

info = zeros(size(x));
% Iterate
while   ~stop  
  if  ng <= tolg,  stop = 1;  
  else
    [h info] = linhuber(-J,r,[],[],[gamma mu Htype],info);
    nh = norm(h);   nx = tolx + norm(x);
    if  nh <= tolx*nx,  stop = 2; end 
  end 
  if  ~stop
    xnew = x + h;   h = xnew - x;   dL = f - info.pf(1);
    [fn Sn rn Jn gn] = huberobj(fun,xnew,gamma,Htype,varargin{:});
    neval = neval + 1;  df = f - fn;  
    if  (dL > 0) & (df > 0)               % Update x and modify mu
      nit = nit + 1;
      x = xnew;   f = fn;  J = Jn;  r = rn;  S = Sn;  g = gn;
      ng = norm(g,inf);
      mu = mu * max(1/3, 1 - (2*df/dL - 1)^3);   nu = 2;
      if  Trace
        X(:,nit) = xnew;   perf(:,nit) = [fn ng mu]'; end
    else                                  % Same  x, increase  mu
      mu = mu*nu;  nu = 2*nu; 
    end
    if  neval == maxeval,  stop = 3; end
  end  
end
%  Set return values
if  Trace
  ii = 1 : nit;  X = X(:,ii);   
  perf = struct('f',perf(1,ii), 'ng',perf(2,ii), 'mu',perf(3,ii));
else,  X = x;  end
info = [f  ng  nh  mu/max(diag(A))  nit-1  neval  stop];