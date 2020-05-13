function  [X,info,perf,B] = smarquardt(fun, x0, opts, B0, varargin)
%SMARQUARDT  Secant version of Levenberg-Marquardt's method for least 
% squares: Find  xm = argmin{f(x)} , where  x = [x_1, ..., x_n]  and
% f(x) = 0.5 * sum(r_i(x)^2) .
% The functions  r_i(x) (i=1,...,m)  must be given by a MATLAB
% function with declaration
%            function  r = fun(x, p1,p2,...)
% p1,p2,... are parameters of the function.  In connection with nonlinear
% data fitting they may be arrays with coordinates of the data points.
%
% Call    [X, info] = smarquardt(fun, x0)
%         [X, info] = smarquardt(fun, x0, opts)
%         [X, info] = smarquardt(fun, x0, opts, B0, p1,p2,...)
%         [X, info, perf] = smarquardt(.....)
%         [X, info, perf, B] = smarquardt(.....)
%
% Input parameters
% fun  :  Handle to the function.
% x0   :  Starting guess for  xm .
% opts :  Either a struct with fields 'tau', 'tolg', 'tolx', 'maxeval'
%         and 'relstep', or a vector with the values of these options,
%         opts = [tau  tolg  tolx  maxeval  relstep].
%         tau  used in starting value for Marquardt parameter: 
%                  mu = tau * max( [B0'*B0]_(ii) ) ,  
%              where B0 is an approximate Jacobian at x0 .
%         relstep  "relative" step length for difference approximations.
%         The other options are used in stopping criteria:
%             ||B(x)'*r(x)||_inf <= tolg                   or 
%             ||dx||_2  <= tolx*(tolx + ||x||_2)           or
%             no. of function evaluations exceeds  maxeval .
%         Default  tau = 1e-3,  tolg = 1e-4,  tolx = 1e-8, 
%             maxeval = 100 + 10*n,  relstep = 1e-7.  
%         If the input opts has less than 5 elements, it is augmented by 
%         the default values.  Also, zeros and negative elements are
%         replaced by the default values.
% B0   :  Initial approximation to J.  If  B0 is not given, a forward 
%         difference approximation to it is used.
% p1,p2,..  are passed dirctly to the function FUN .    
%
% Output parameters
% X    :  If  perf  is present, then array, holding the iterates
%         columnwise.  Otherwise, computed solution vector.
% info :  Performance information, vector with 7 elements:
%         info(1:4) = final values of 
%           [f(x)  ||B(x)'*r(x)||inf  ||dx||2  mu/max( [B(x)'*B(x)]_(ii) )] .
%         info(5:6) = no. of iteration steps and function evaluations
%         info(6) = 1 :  Stopped by small gradient
%                   2 :  Stopped by small x-step
%                   3 :  No. of iteration steps exceeded 
%                  -1 :  x is not a real valued vector
%                  -2 :  r is not a real valued column vector
%                  -4 :  Dimension mismatch in x, r, B0
%                  -5 :  Overflow during computation
%                  -6 :  Error in approximate Jacobian
%         info(7) = no. of iterations.
% perf :  Struct with fields
%          f :  values of  f(xk) ,
%         ng :  values of  || B'*r(x) ||_inf ,
%         mu :  values of Marquardt damping parameter.
% B    :  Computed approximation to Jacobian at the solution.

% Version 10.11.08.  hbn(a)imm.dtu.dk

% Check parameters and function call
f = NaN;  ng = NaN;  perf = [];  B = [];
info = zeros(1,7);
if  nargin < 2,  stop = -1;
else
  [stop x n] = checkx(x0);   
  if  ~stop
    [stop f r] = checkrJ(fun,x0,varargin{:});  info(6) = 1;
    if  ~stop
      %  Finish initialization
      if  nargin < 3,  opts = []; end
      opts  = checkopts('smarquardt', opts);  % use default options where required
      tau = opts(1);    tolg = opts(2);  tolx = opts(3);  relstep = opts(5);  
      if  opts(4) > 0,  maxeval = opts(4); else,  maxeval = 100 + 10*n; end
      % Jacobian
      if  nargin > 3  % B0 is given
        sB = size(B0);  m = length(r);
        if  sum(sB) == 0  % placeholder
          [stop B] = Dapprox(fun,x,relstep,r,varargin{:});  
          info(6) = info(6) + n;
        elseif  any(sB ~= [m n]),  stop = -4;
        else,                      B = B0;   end
      else
        [stop B] = Dapprox(fun,x,relstep,r,varargin{:});  
        info(6) = info(6) + n;
      end
      % Check gradient and J'*J
      if  ~stop
        g = B'*r;   ng = norm(g,inf);  A = B'*B;
        if  isinf(ng) | isinf(norm(A(:),inf)),  stop = -5; end 
      end
    end
  end
end
if  stop
  X = x0;  info([1:5 7]) = [f ng  0  tau  0 stop];
  return
end

% Finish initialization
mu = tau * max(diag(A));    
Trace = nargout > 2;
if  Trace
  o = ones(1, maxeval);  
  X = x * o;  perf = [f; ng; mu] * o;
end 

% Iterate
k = 1;   nu = 2;   nh = 0;
ng0 = ng;
ku = 0;   % direction of last update
kit = 0;  % no. of iteration steps

while  ~stop
  if  ng <= opts(2),  stop = 1; 
  else
    [h mu] = geth(A,g,mu);
    nh = norm(h);   nx = tolx + norm(x);
    if  nh <= tolx*nx,  stop = 2; end 
  end 
  if  ~stop
    xnew = x + h;    h = xnew - x;  
    [stop fn rn] = checkrJ(fun,xnew,varargin{:});  info(6) = info(6)+1;
    if  ~stop
      % Update  B
      ku = mod(ku,n) + 1; 
      if  abs(h(ku)) < .8*norm(h)  % extra step
        xu = x;
        if  x(ku) == 0,  xu(ku) = opts(5)^2;
        else,            xu(ku) = x(ku) + opts(5)*abs(x(ku)); end
        [stop fu ru] = checkrJ(fun,xu,varargin{:});  info(6) = info(6)+1;
        if  ~stop
          hu = xu - x;
          B = B + ((ru - r - B*hu)/norm(hu)^2) * hu';
        end
      end
      B = B + ((rn - r - B*h)/norm(h)^2) * h'; 
      k = k + 1;
      dL = (h'*(mu*h - g))/2;   
      if  length(rn) ~= length(r)
        df = f - fn;
      else  % more accurate
        df = ( (r - rn)' * (r + rn) )/2; 
      end 
      if  (dL > 0) & (df > 0)               % Update x and modify mu
        kit = kit + 1;   
        x = xnew;   f = fn;  r = rn;
        mu = mu * max(1/3, 1 - (2*df/dL - 1)^3);   nu = 2;
        if  Trace
          X(:,kit+1) = x;   perf(:,kit+1) = [fn norm(B'*rn,inf) mu]'; end
      else  % Same  x, increase  mu
        mu = mu*nu;  nu = 2*nu; 
      end 
      if  info(5) > maxeval,  stop = 3; 
      else
        g = B'*r;  ng = norm(g,inf);      A = B'*B;
        if  isinf(ng) | isinf(norm(A(:),inf)),  stop = -5; end
      end
    end  
  end
end
%  Set return values
if  Trace
  ii = 1 : kit+1;  X = X(:,ii);   
  perf = struct('f',perf(1,ii), 'ng',perf(2,ii), 'mu',perf(3,ii));
else,  X = x;  end
if  stop < 0,  tau = NaN;  else,  tau = mu/max(diag(A)); end
info([1:5 7]) = [f  ng  nh  tau  kit stop];