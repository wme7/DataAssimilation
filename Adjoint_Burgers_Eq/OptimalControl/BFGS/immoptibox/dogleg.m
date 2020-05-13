function  [X, info, perf] = dogleg(fun, x0, opts, varargin)
%DOGLEG Powell's dog-leg method for least squares.
% Find  xm = argmin{f(x)} , where  x  is an n-vector and
%     f(x) = 0.5 * sum(r_i(x)^2) .
% The functions  r_i(x) (i=1,...,m) and the Jacobian matrix  J(x)    
% (with elements  J(i,j) = Dr_i/Dx_j ) must be given by a MATLAB
% function with declaration
%            function  [r, J] = fun(x, p1,p2,...)
% The gradient of  f  is  g = J' * r.
% p1,p2,... are parameters of the function.  In connection with nonlinear 
% data fitting they may be arrays with coordinates of the data points.
%  
% Call    [X, info] = dogleg(fun, x0)
%         [X, info] = dogleg(fun, x0, opts, p1,p2,...)
%         [X, info, perf] = dogleg(.....)
%
% Input parameters
% fun  :  Handle to the function.
% x0   :  Starting guess for  xm .
% opts :  Either a struct with fields 'Delta', 'tolg', 'tolx', 'tolr' and 
%         'maxeval', or a vector with the values of these options,
%         opts = [Delta  tolg  tolx  tolr  maxeval].
%         Delta    initial trust region radius.
%         The other options are used in stopping criteria:
%             ||g(x)||_inf <= tolg                         or 
%             ||dx||_2  <= tolx*(tolx + ||x||_2)           or
%             ||r(x)||_inf <= tolr                         or
%             no. of function evaluations exceeds  maxeval .
%         Default  Delta = [0.1(1+||x0||),  tolg = 1e-4,  tolx = 1e-8,  
%             tolr = 1e-6,  maxeval = 100. 
%         If the input opts has less than 5 elements, it is augmented by 
%         the default values.  Also, zeros and negative elements are
%         replaced by the default values.
% p1,p2,..  are passed dirctly to the function FUN .    
%
% Output parameters
% X    :  If  perf  is present, then array, holding the iterates  xk
%         columnwise.  Otherwise, computed solution vector.
% info :  Performance information, vector with 7 elements:
%         info(1:4) = final values of  [f(x)  ||g||_inf  ||dx||_2  Delta] 
%         info(5:6) = no. of iteration steps and function evaluations
%         info(7) = 1 :  Stopped by small gradient
%                   2 :  Stopped by small x-step
%                   3 :  Stopped by small r-vector
%                   4 :  No. of iteration steps exceeded
%                  -1 :  x is not a real valued vector
%                  -2 :  r is not a real valued column vector
%                  -3 :  J is not a real valued matrix
%                  -4 :  Dimension mismatch in x, r, J
%                  -5 :  Overflow during computation 
% perf :  Struct with fields
%             f :  values of  f(xk) ,
%            ng :  values of  || g(xk) ||_inf ,
%         Delta :  values of the trust region radius.

% Version 10.10.08.  hbn(a)imm.dtu.dk

% Check parameters and function call
if  nargin < 2,  stop = -1;
else
  [stop x n] = checkx(x0);   
  if  ~stop,  [stop f r J] = checkrJ(fun,x0,varargin{:});  k = 1; end
end
if  ~stop
  g = -(J'*r);   ngi = norm(g,inf);  % negative gradient and its norm
  if  isinf(ngi),  stop = -5; end
else
  f = NaN;  ngi = NaN;
end
if  stop
  X = x0;  perf = [];  info = [f  ngi  0  opts(1)  0 1  stop];
  return
end
  
% Finish initialization
if  nargin < 3,  opts = []; end
opts  = checkopts('dogleg', opts);  % use default options where required
if  opts(1) == 0,  opts(1) = 0.1 * (1 + norm(x)); end
Delta = opts(1);    tolg = opts(2);  tolx = opts(3);  tolr = opts(4); 
maxeval = opts(5);
Trace = nargout > 2;
if  Trace
  o = ones(1, maxeval);  
  X = x * o;  perf = [f; ngi; Delta] * o;
end 

nu = 2;   nx = norm(x);  fact = 1;  stop = 0;
reduce = 0;   nstep = 0;  kit = 0;

% Iterate
while   ~stop 
  if      isinf(ngi),                  stop = -5;
  elseif  ngi <= tolg,                 stop = 1;  
  elseif  Delta <= tolx*(tolx + nx),   stop = 2;  
  elseif  norm(r,inf) <= tolr,         stop = 3; 
  elseif  k >= maxeval,                stop = 4;
  else
    if  fact  % Factorize and compute hGN
      [Q R] = qr(J,0);  [U S V] = svd(R);
      s = diag(S);  i = find(s > 100*eps*s(1));
      Qr = -(r'*Q)';   UQr = U(:,i)'*Qr;
      hGN = V(:,i)*(UQr ./ s(i));  
      nhGN = norm(hGN);  fact = 0;
    end
    if  nhGN > Delta  % Include gradient
      nh = Delta;
      ng = norm(g);  alpha = (ng / norm(J*g))^2;
      gn = alpha*g;  ngn = alpha*ng;
      if  ngn >= Delta
        h = (Delta/ngn) * gn;
        dLpre = Delta*(2*ngn - Delta)/(2*alpha);
      else  %  Dog leg
        b = hGN - gn;  bb = b'*b;  gb = gn'*b;
        c = (Delta + ngn)*(Delta - ngn);
        if  gb > 0
          beta = c / (gb + sqrt(gb^2 + c * bb));
        else
          beta = (sqrt(gb^2 + c * bb) - gb)/bb;
        end
        h = gn + beta*b;
        dLpre = .5*alpha*(1 - beta)^2*ng^2 + beta*(2-beta)*f;
      end
    else
      h = hGN;  nh = nhGN;
      dLpre = f;  
      if  nh <= tolx*(tolx + nx),  stop = 2; end
    end
  end
  if  ~stop
    xnew = x + h;   h = xnew - x;   nstep = norm(h);
    dL = f - .5*norm(r + J*h)^2;
    [stop fn rn Jn] = checkrJ(fun,xnew,varargin{:});  k = k + 1;
    if  ~stop
      df = f - fn;
      if  df > 0 & dL > 0                        % Update x 
        kit = kit + 1;  updx = 1;
        x = xnew;   f = fn;  J = Jn;  r = rn;   fact = 1;
        g = -(J'*r);   ngi = norm(g,inf);
        rho = dL / df; 
      else,  rho = -1;  updx = 0; end
      % Update Delta
      if  abs(rho-1) < 0.2 & nh > Delta/3 & reduce <= 0
        Delta = 3*Delta;   nu = 2;  reduce = 0;
      elseif  rho < 0.25
        Delta = Delta/nu;
        nu = 2*nu;  reduce = 2;
      else
        reduce = reduce - 1;
      end       
      if  Trace && updx  % store iterate
        X(:,kit+1) = xnew;   perf(:,kit+1) = [fn ngi Delta]'; 
      end
      if  k > maxeval,  stop = 3; end 
    end
  end   
end
%  Set return values
if  Trace
  ii = 1 : kit+1;  X = X(:,ii);   
  perf = struct('f',perf(1,ii), 'ng',perf(2,ii), 'Delta',perf(3,ii));
else,  X = x;  end
if  stop < 0,  f = NaN;  ngi = NaN; end
info = [f  ngi  nstep  Delta  kit k  stop];