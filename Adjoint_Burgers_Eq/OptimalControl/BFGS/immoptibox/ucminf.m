function  [X, info, perf, D] = ucminf(fun, x0, opts, D0, varargin)
%UCMINF  BFGS method for unconstrained nonlinear optimization:
% Find  xm = argmin{f(x)} , where  x  is an n-vector and the scalar
% function  f  with gradient  g  (with elements  g(i) = Df/Dx_i )
% must be given by a MATLAB function with with declaration
%            function  [f, g] = fun(x,p1,p2,...)
% p1,p2,... are parameters of the function.
% 
% Call  [X, info] = ucminf(fun, x0)
%       [X, info] = ucminf(fun, x0, opts)
%       [X, info] = ucminf(fun, x0, opts, D0, p1,p2,...)
%       [X, info, perf] = ucminf(.....)
%       [X, info, perf, D] = ucminf(.....)
%
% Input parameters
% fun  :  Handle to the function.
% x0   :  Starting guess for  xm .
% opts :  Either a struct with fields 'Delta', 'tolg', 'tolx' and 'maxeval'
%         or a vector with the values of these options,
%         opts = [Delta  tolg  tolx  maxeval].
%         Delta :  Expected length of initial step.
%         The other options are used in stopping criteria:
%             ||g(x)||_inf <= tolg                     or 
%             ||dx||2 <= tolx*(tolx + ||x||_2)         or
%             no. of function evaluations exceeds  maxeval .
%         Default  Delta = 1,  tolg = 1e-4, tolx = 1e-8,  maxeval = 100.
%         If the input opts has less than 4 elements, it is augmented by 
%         the default values.  Also, zeros and negative elements are
%         replaced by the default values.
% D0   :  If present, then approximate inverse Hessian at  x0 .
%         Otherwise, D0 := I.
% p1,p2,..  are passed dirctly to the function FUN .
%
% Output parameters
% X    :  If  perf  is present, then array, holding the iterates
%         columnwise.  Otherwise, computed solution vector.
% info :  Performance information, vector with 6 elements:
%         info(1:3) = final values of [f(x)  ||g||_inf  ||dx||_2] 
%         info(4:5) = no. of iteration steps and evaluations of (f,g)
%         info(6) = 0 :  Line search failed
%                   1 :  Stopped by small gradient
%                   2 :  Stopped by small x-step
%                   3 :  No. of iteration steps exceeded 
%                   4 :  Stopped by zero step.
%                  -1 :  x is not a real valued vector
%                  -2 :  f is not a real valued scalar
%                  -3 :  g is not a real valued vector 
%                  -4 :  x and g have different lengths
%                  -6 :  D0 is not a symmetric, pos. def. n*n-matrix% 
% perf :  Struct with fields
%             f :  values of  f(xk) ,
%            ng :  values of  || g(xk) ||_inf ,
%         Delta :  values of trust region radius.
%            am :  am from line search,
%         slope :  P'(am) from line search,
%         neval :  no. of fct. evals in current line search.
% D   :   Array holding the approximate inverse Hessian at the computed 
%         minimizer.
%
% Line search is performed by the immoptibox function LINESEARCH.

% Version 10.09.29.  hbn(a)imm.dtu.dk

% Initial check
if  nargin < 2,  error('Too few input parameters'), end

% Check OPTS
if  nargin < 3,  opts = []; end
opts  = checkopts('ucminf', opts);  % use default options where required
Delta = opts(1);  tolg = opts(2);  tolx = opts(3);  maxeval = opts(4);

% Check parameters and function call
[stop x n] = checkx(x0);   
if  ~stop
  [stop f g] = checkfgH(fun,x0,varargin{:}); 
  if  ~stop  % Initial inverse Hessian
    if  nargin > 3 & ~isempty(D0)
      [stop D] = checkD(n,D0);  fst = 0;
    else,   D = eye(n);  fst = 1;    end
  end
else,  f = NaN;  end

if  stop
  X = x0;  perf = [];  D = [];  info = [f(1)  NaN  0  0  1  stop];
  return
end

%  Finish initialization
k = 1;   neval = 1;   ng = norm(g,inf);
Trace = nargout > 2;
if  Trace
  o = ones(1, maxeval);  
  X = x * o;  perf = [f; ng; Delta; 0; 0; 0] * o;
end
if  ng <= tolg,  stop = 1;  nh = 0;
else
  h = zeros(size(x));  nh = 0;
  ngs = repmat(ng,1,3);
  lsopts = [1 .05 .99 5 2];  % LINESEARCH options
end

more = 1;
while  ~stop & more
  %  Previous values
  xp = x;   gp = g;   fp = f;   nx = norm(x);
  ngs = [ngs(2:3) ng];
  h = D*(-g(:));   nh = norm(h);   red = 0; 
  if  nh <= tolx*(tolx + nx),  stop = 2;  
  else
    if  fst | nh > Delta  % Scale to ||h|| = Delta
      h = (Delta / nh) * h;   nh = Delta;   
      fst = 0;  red = 1;
    end
    %  Line search
    [x f g  linfo] = linesearch(fun,x,f,g, h, lsopts, varargin{:}); 
    neval = neval + linfo(3);  
    if  linfo(1) <= 0
      stop = linfo(1);  more = 0;  % something wrong
    else
      k = k+1;
      if  linfo(1) < 1  % Reduce Delta
        Delta = .35 * Delta;
      elseif   red & (linfo(2) > .7)  % Increase Delta
        Delta = 3*Delta;      
      end 
      %  Update ||g||
      ng = norm(g,inf);
      if  Trace
        X(:,k) = x(:); 
        perf(:,k) = [f; ng; Delta; linfo(1); dot(h,g); linfo(3)]; 
      end
      h = x - xp;   nh = norm(h);
      if  nh == 0,
        stop = 4; 
      else
        y = g - gp;   yh = dot(y,h);
        if  yh > sqrt(eps) * nh * norm(y)
          %  Update  D
          v = D*y(:);   yv = dot(y,v);
          a = (1 + yv/yh)/yh;   w = (a/2)*h(:) - v/yh;
          D = D + w*h' + h*w';
        end  % update D
        %  Check stopping criteria
        thrx = tolx*(tolx + norm(x));
        if      ng <= tolg,                 stop = 1;
        elseif  nh <= thrx,                 stop = 2;
        elseif  neval >= maxeval,           stop = 3; 
        else,   Delta = max(Delta, 2*thrx);  end
      end  
    end  % Nonzero h
  end % nofail
end  % iteration

%  Set return values
if  Trace
  ii = 1 : k;  X = X(:,ii);   
  perf = struct('f',perf(1,ii), 'ng',perf(2,ii), 'Delta',perf(3,ii), ...
    'am',perf(4,ii), 'slope',perf(5,ii), 'neval',perf(6,ii));
else,  X = x;  end
info = [f  ng  nh  k-1  neval  stop];

% ==========  auxiliary function  =================================

function  [err, D] = checkD(n,D0)
% Check given inverse Hessian
D = D0;   sD = size(D);  err = 0;
if  any(sD - n),  err = -6;  return, end
% Check symmetry
dD = D - D';   ndD = norm(dD(:),inf);
if  ndD > 10*eps*norm(D(:),inf),  err = -6;  return, end
if  ndD,  D = (D + D')/2; end  % Symmetrize      
[R p] = chol(D);
if  p,  err = -6;  end