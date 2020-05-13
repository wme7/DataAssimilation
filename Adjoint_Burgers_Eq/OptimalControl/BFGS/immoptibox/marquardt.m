function  [X, info, perf] = marquardt(fun, x0, opts, varargin)
%MARQUARDT  Levenberg-Marquardt's method for least squares.
% Find  xm = argmin{f(x)} , where  x  is an n-vector and
% f(x) = 0.5 * sum(r_i(x)^2) .
% The functions  r_i(x) (i=1,...,m) and the Jacobian matrix  J(x)  
% (with elements  J(i,j) = Dr_i/Dx_j ) must be given by a MATLAB
% function with declaration
%            function  [r, J] = fun(x,p1,p2,...)
% p1,p2,... are parameters of the function.  In connection with 
% nonlinear data fitting they may be arrays with coordinates of 
% the data points.
%  
% Call    [X, info] = marquardt(fun, x0)
%         [X, info] = marquardt(fun, x0, opts, p1,p2,...)
%         [X, info, perf] = marquardt(.....)
%
% Input parameters
% fun  :  Handle to the function.
% x0   :  Starting guess for  xm .
% opts :  Either a struct with fields 'tau', 'tolg', 'tolx' and 'maxeval', 
%         or a vector with the values of these options,
%         opts = [tau  tolg  tolx  maxeval].
%         tau    used in starting value for Marquardt parameter: 
%             mu = tau * max( [J(x0)'*J(x0)]_(ii) )  .
%         The other options are used in stopping criteria:
%             ||g||_inf <= tolg                            or 
%             ||dx||_2  <= tolx*(tolx + ||x||_2)           or
%             no. of function evaluations exceeds  maxeval .
%         Default  tau = 1e-3,  tolg = 1e-4,  tolx = 1e-8, maxeval = 100. 
%         If the input opts has less than 4 elements, it is augmented by 
%         the default values.  Also, zeros and negative elements are
%         replaced by the default values.
% p1,p2,..  are passed directly to the function FUN . 
%
% Output parameters
% X    :  If  perf  is present, then array, holding the iterates
%         columnwise.  Otherwise, computed solution vector.
% info :  Performance information, vector with 7 elements:
%         info(1:4) = final values of 
%             [f(x)  ||g(x)||inf  ||dx||2  mu/max( [J(x)'*J(x)]_(ii) )] .
%         info(5:6) = no. of iteration steps and function evaluations.
%         info(7) = 1 :  Stopped by small gradient
%                   2 :  Stopped by small x-step
%                   3 :  No. of function evaluations exceeded 
%                  -1 :  x is not a real valued vector
%                  -2 :  r is not a real valued column vector
%                  -3 :  J is not a real valued matrix
%                  -4 :  Dimension mismatch in x, r, J
%                  -5 :  Overflow during computation 
% perf :  Struct with fields
%          f :  values of  f(xk) ,
%         ng :  values of  || g(xk) ||_inf ,
%         mu :  values of Marquardt damping parameter.

% Version 10.11.03.  hbn(a)imm.dtu.dk

% Check parameters and function call
if  nargin < 2,  stop = -1;
else
  [stop x n] = checkx(x0);   
  if  ~stop,  [stop f r J] = checkrJ(fun,x0,varargin{:});  k = 1; end
end
if  ~stop
  g = J'*r;   ng = norm(g,inf);  A = J'*J;
  if  isinf(ng) | isinf(norm(A(:),inf)),  stop = -5; end
else
  f = NaN;  ng = NaN;
end
if  stop
  X = x0;  perf = [];  info = [f  ng  0  opts(1)  0  1 stop];
  return
end

if  nargin < 3,  opts = []; end
opts  = checkopts('marquardt', opts);  % use default options where required
tau = opts(1);    tolg = opts(2);  tolx = opts(3);  maxeval = opts(4); 

%  Finish initialization
mu = tau * max(diag(A));
Trace = nargout > 2;
if  Trace
  o = ones(1, maxeval);  
  X = x * o;  perf = [f; ng; mu] * o;
end 
nu = 2;   nh = 0;   stop = 0;  kit = 0;

% Iterate
while   ~stop  
  if  ng <= tolg,  stop = 1;
  else
    [h mu] = geth(A,g,mu); 
    nh = norm(h);   nx = tolx + norm(x);
    if  nh <= tolx*nx,  stop = 2; end 
  end 
  if  ~stop
    xnew = x + h;   h = xnew - x;   dL = (h'*(mu*h - g))/2; 
    [stop fn rn Jn] = checkrJ(fun, xnew, varargin{:});  k = k + 1;
    if  ~stop 
      if  length(rn) ~= length(r)
        df = f - fn;
      else  % more accurate
        df = ( (r - rn)' * (r + rn) )/2; 
      end  
      if  (dL > 0) & (df > 0)               % Update x and modify mu
        kit = kit + 1;
        x = xnew;   f = fn;  J = Jn;  r = rn; 
        A = J'*J;   g = J'*r;   ng = norm(g,inf);
        mu = mu * max(1/3, 1 - (2*df/dL - 1)^3);   nu = 2;
      if  Trace
        X(:,kit+1) = xnew;   perf(:,kit+1) = [fn ng mu]'; end
      else                                  % Same  x, increase  mu
        mu = mu*nu;  nu = 2*nu; 
      end
      if      k > maxeval,                        stop = 3; 
      elseif  isinf(ng) | isinf(norm(A(:),inf)),  stop = -5; end
    end
  end   
end
%  Set return values
if  Trace
  ii = 1 : kit+1;  X = X(:,ii);   
  perf = struct('f',perf(1,ii), 'ng',perf(2,ii), 'mu',perf(3,ii));
else,  X = x;  end
if  stop < 0,  f = NaN;  ng = NaN; end
info = [f  ng  nh  mu/max(diag(A))  kit k  stop];