function  [X, info, perf] = dampnewton(fun, x0, opts, varargin)
%DAMPNEWTON  Levenberg-Marquardt damped Newton method for unconstrained 
% optimization:  Find  xm = argmin{f(x)} , where  x  is an n-vector and the
% scalar function  f  with gradient  g  (with elements  g_i = Df/Dx_i )
% and Hessian  H  (with elements  H_{j,k} = D^2f/(Dx_j Dx_k) ) must be 
% given by a Matlab function with with declaration
%            function  [f, g, H] = fun(x, p1,p2,...)
% p1,p2,... are parameters of the function.
%  
% Call    [X, info] = dampnewton(fun, x0)
%         [X, info] = dampnewton(fun, x0, opts, p1,p2,...)
%         [X, info, perf] = dampnewton(.....)
%
% Input parameters
% fun  :  Handle to the function.
% x0   :  Starting guess for  xm .
% opts :  Either a struct with fields 'tau', 'tolg', 'tolx' and 'maxeval'
%         or a vector with the values of these options,
%         opts = [tau  tolg  tolx  maxeval].
%         tau  is used in starting value for damping parameter: 
%             mu = tau) * ||H0||_inf, where H0 is the  Hessian of f at x0.
%         The other options are used in stopping criteria:
%             ||g(x)||_inf <= tolg                     or 
%             ||dx||2 <= tolx*(tolx + ||x||_2)         or
%             no. of function evaluations exceeds  maxeval .
%         Default  tau = 1e-3,  tolg = 1e-4,  tolx = 1e-8,  maxeval = 100.
%         If the input opts has less than 4 elements, it is augmented by 
%         the default values.  Also, zeros and negative elements are
%         replaced by the default values.
% p1,p2,..  are passed directly to the function FUN .    
%
% Output parameters
% X    :  If  perf  is present, then array, holding the iterates  xk
%         columnwise.  Otherwise, computed solution vector.
% info :  Performance information, vector with 7 elements:
%         info(1:4) = final values of 
%             [f(x)  ||g||_inf  ||dx||_2  mu/||H(x)||_inf] .
%         info(5:6) = no. of iteration steps and function evaluations.
%         info(7) = 1 : Stopped by small gradient
%                   2 :  Stopped by small x-step
%                   3 :  No. of function evaluations exceeded 
%                  -1 :  x is not a real valued vector
%                  -2 :  f is not a real valued scalar
%                  -3 :  g is not a real valued vector or
%                        H is not a real valued matrix
%                  -4 :  Dimension mismatch in x, g, H
%                  -5 :  H is not symmetric
% perf :  Struct with fields
%            f :  values of  f(xk) ,
%           ng :  values of  || g(xk) ||_inf ,
%           mu :  values of the damping  mu(xk) .

% Version 10.10.08.  hbn(a)imm.dtu.dk

% Check parameters and function call
if  nargin < 2,  stop = -1;
else
  [stop x n] = checkx(x0);   
  if  ~stop
    [stop f g H] = checkfgH(fun,x0,varargin{:});
  end
end
if  stop
  X = x0;  perf = [];  info = [repmat(NaN,1,4) 0  1  stop];
  return
end
if  nargin < 3,  opts = []; end
opts  = checkopts('dampnewton', opts);  % use default options where required
tau = opts(1);  tolg = opts(2);  tolx = opts(3);  maxeval = opts(4);

%  Finish initialization
ng = norm(g,inf);   mu = tau * norm(H,inf);  
if  mu == 0  % zero initial Hessian.  Steepest descent direction
  mu = ng / max(norm(x), sqrt(eps));
end

Trace = nargout > 2;
if  Trace
  o = ones(1, maxeval);  
  X = x * o;  perf = [f; ng; mu] * o;
end 
nu = 2;   nh = 0;   stop = 0;

info(5:7) = [0 1 0];  kit = 0;
% Iterate
while   ~stop 
  if  ng <= tolg,  stop = 1;  
  else
    [h mu] = geth(H,g,mu);
    nh = norm(h);   nx = tolx + norm(x);
    if  nh <= tolx*nx,  stop = 2; end 
  end 
  if  ~stop
    xnew = x + h;   h = xnew - x;   
    [stop fn gn Hn] = checkfgH(fun,xnew,varargin{:});  
    info(6) = info(6) + 1;
    if  ~stop
      df = f - fn;   accept = 0;
      if  df > 0  
        accept = 1;  dL = (h'*(mu*h - g))/2; 
      elseif  fn <= f + abs(f)*(1 + 100*eps)  % Try gradient
        df = (g + gn)'*(g - gn);
        if  df > 0,  accept = 2; end
      end        
      if  accept                     % Update x and modify mu
        kit = kit + 1;
        if  accept == 1 & dL > 0
          mu = mu * max(1/3, 1 - (2*df/dL - 1)^3);  nu = 2;
        else
          mu = mu*nu;  nu = 2*nu;  
        end
        x = xnew;   f = fn;  g = gn;  H = Hn;   ng = norm(g,inf);
        if  Trace
          jtr = kit + 1;
          X(:,jtr) = xnew;   perf(:,jtr) = [f ng mu]'; end
      else                                  % Same  x, increase  mu
        mu = mu*nu;  nu = 2*nu; 
      end
      if  info(6) == maxeval, stop = 3;  end
    end
  end   
end
%  Set return values
if  Trace
  jj = 1:jtr;  X = X(:,jj);   
  perf = struct('f',perf(1,jj), 'ng',perf(2,jj), 'mu',perf(3,jj));
else,  X = x;  end
if  stop < 0,  f = NaN;  ng = NaN; end
info = [f  ng  nh  mu / norm(H,inf)  kit  info(6)  stop];