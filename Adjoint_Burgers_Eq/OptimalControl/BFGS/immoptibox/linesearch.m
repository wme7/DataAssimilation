function  [xn,fn,gn,info,perf] = ...
                   linesearch(fun, x,f,g, h, opts, varargin)
%LINESEARCH  Find  am = argmin_{a > 0}{ P(a) = f(x+a*h) } , where  x  and  
% h  are given n-vectors and the scalar function  f  and its gradient  g  
% (with elements  g(i) = Df/Dx_i ) must be given by a MATLAB function with 
% declaration
%            function  [f, g] = fun(x,p1,p2,...)
% p1,p2,... are parameters of the function.
%
%  Call  [xn,fn,gn,info] = linesearch(fun,x,f,g,h)
%        [xn,fn,gn,info] = linesearch(fun,x,f,g,h,opts,p1,p2,...)
%        [xn,fn,gn,info,perf] = linesearch(......)
%
% Input parameters
% fun  :  Handle to the function.
% x    :  Current x.
% f,g  :  f(x) and g(x).
% h    :  Step vector.
% opts :  Either a struct with fields 'choice', 'cp1', 'cp2', 'maxeval',
%         'amax'  or a vector with the values of these options,
%         opts = [choice  cp1  cp2  maxeval  amax].
%         choice = 0 :  exact line search.  
%                       Otherwise: soft line search  (Default).
%         cp1, cp2   :  options for stopping criteria.
%              choice = 0:  |P'(a)| <= cp1*|P'(0)|  or  c-b <= cp2*c , 
%                           where  [b,c] is the current interval for a.  
%                           Default  cp1 = cp2 = 1e-3 .
%              Otherwise:  P(a) <= P(0) + a*cp1*P'(0)  and  
%                          P'(a) >= cp2*P'(0).  
%                          Default  cp1 = 1e-3,  cp2 = 0.99 .
%         maxeval    :  Maximum number of function evaluations.
%                       Default: maxeval = 10.
%         amax       :  Maximal allowable step.   Default  amax = 10.
% p1,p2,..  are passed dirctly to the function FUN .    
%
% Output parameters
% xn    :  x + am*h
% fn,gn :  f(xn) and g(xn).
% info  :  Performance information, vector with 3 elements
%          info(1) >  0 : am.  Successfull call
%                  =  0 : h is not downhill or it is so large and maxeval
%                         so small, that a better point was not found.
%                  = -1 : x is not a real valued vector.
%                  = -2 : f is not a real valued scalar.
%                  = -3 : g or h is not a real valued vector.
%                  = -4 : g or h has different length from x.
%          info(2) = slope ratio at the solution,  P'(am)/P'(0) .
%          info(3) = number of function evaluations used.
% perf  :  Struct with fields
%          alpha :  values of  a ,
%            phi :  values of  P(a) = f(x+a*h) ,
%          slope :  values of  P'(a) .

% Version 10.09.28.  hbn(a)imm.dtu.dk

% Initial check
if  nargin < 5,  error('Too few input parameters'), end

% Check OPTS
if  nargin < 6 | isempty(opts),  opts = 1; end
opts  = checkopts('linesearch', opts);  % use default options where required
% if  opts(1) == 0,  opts = checkopts(opts, [0  1e-3  1e-3  10 10]); 
% else,              opts = checkopts(opts, [1  1e-3  0.99  10 10]); end
choice = opts(1) ~= 0;
cp1 = opts(2);  cp2 = opts(3);  maxeval = opts(4);  amax = opts(5);

% Default return values and simple checks
xn = x;  fn = f;   gn = g;  info = [0 1 0];
[info(1)  n] = check(x,f,g,h);
if  info(1),  return,  else,  stop = 0; end

x = x(:);   h = h(:);  % both are treated as column vectors
% Check descent condition
f0 = fn;  df0 = dot(h,gn); 
if  df0 >= -10*eps*norm(h)*norm(gn)  % not significantly downhill
  info(1) = 0;  return
end

Trace = nargout > 4;
if  Trace
  o = ones(1, maxeval);  
  X = x * o;  perf = [0; f0; df0] * o;
end 

% Finish initialization 
if  choice  % soft line search
  slope0 = cp1*df0;   slopethr = cp2*df0;
else  % exact line search
  slope0 = 0;   slopethr = cp1*abs(df0);
end

% Get an initial interval for am
a = 0;  fa = fn;  dfa = df0;  stop = 0;
b = min(1, amax);
while  ~stop
  [stop fb g] = checkfgH(fun,x+b*h,varargin{:});   info(3) = info(3)+1;
  if  stop,  info(1) = stop;
  else
    dfb = dot(g,h);
    if  Trace,  perf(:,info(3)) = [b; fb; dfb]; end
    if  fb < f0 + slope0*b  % new lower bound
      info(1:2) = [b dfb/df0];
      if  choice,  a = b;  fa = fb;  dfa = dfb; end
      xn = x + b*h;  fn = fb;  gn = g;
      if  (dfb < min(slopethr,0)) && (info(3) < maxeval) && (b < amax)
        % Augment right hand end
        if  ~choice,  a = b;  fa = fb;  dfa = dfb; end
        if  2.5*b >= amax,  b = amax;  else,  b = 2*b; end
      else,  stop = 1; end
    else,  stop = 1; end
  end
end % phase 1: expand interval

if  stop >= 0  % OK so far.  Check stopping criteria
  stop = (info(3) >= maxeval) | (b >= amax & dfb < slopethr)...  % Cannot improve
         | (choice & (a > 0 & dfb >= slopethr));  % OK
end
if  stop
  if  Trace
    ii = 1 : info(3);
    perf = struct('alpha',perf(1,1:k), 'phi',perf(2,1:k), 'slope',perf(3,1:k)); 
  end
  return
end

% Refine interval.  Use auxiliary array  xfd
xfd = [a b b; fa fb fb; dfa dfb dfb];
while  ~stop
  c = interpolate(xfd,n);
  [stop fc g] = checkfgH(fun,x+c*h,varargin{:});   info(3) = info(3)+1;
  if  stop,  info(1) = stop;
  else
    xfd(:,3) = [c; fc; dot(g,h)];
    if  Trace,  perf(:,info(3)) = xfd(:,3); end
    if  choice  % soft line search
      if  fc < f0 + slope0*c  % new lower bound
        info(1:2) = [c xfd(3,3)/df0];
        xn = x + c*h;  fn = fc;  gn = g;  
        xfd(:,1) = xfd(:,3);
        stop = xfd(3,3) > slopethr;
      else  % new upper bound
        xfd(:,2) = xfd(:,3);  
      end
    else  % exact line search
      if  fc < fn  % better approximant
        info(1:2) = [c xfd(3,3)/df0];
        xn = x + c*h;  fn = fc;  gn = g;  
      end
      if  xfd(3,3) < 0,  xfd(:,1) = xfd(:,3);  % new lower bound
      else,         xfd(:,2) = xfd(:,3);  end  % new upper bound
      stop = abs(xfd(3,3)) <= slopethr...
        | diff(xfd(1,1:2)) < cp2*xfd(1,2);
    end
  end
  stop = stop | info(3) >= maxeval;
end % refine

% Return values
if  Trace
  ii = 1 : info(3);
  perf = struct('alpha',perf(1,ii), 'phi',perf(2,ii), 'slope',perf(3,ii));
end

%============  Auxiliary functions  ========================

function  t = interpolate(xfd,n);
% Minimizer of parabola given by  xfd = [a b; f(a) f(b); f'(a) dummy]
a = xfd(1,1);   b = xfd(1,2);   d = b - a;   df = xfd(3,1);
C = diff(xfd(2,1:2)) - d*df;
if C >= 5*n*eps*b    % Minimizer exists
  A = a - .5*df*(d^2/C);  d = 0.1*d;   
  t = min(max(a+d, A), b-d);  % Ensure significant resuction
else
  t = (a+b)/2;
end

function  [err, n] = check(x,f,g,h)
% Check  x
err = 0;  sx = size(x);   n = max(sx);
if  (min(sx) ~= 1) | ~isreal(x) | any(isnan(x(:))) | isinf(norm(x(:))) 
  err = -1; 
else
  % Check  f
  sf = size(f); 
  if  any(sf ~= 1) | ~isreal(f) | any(isnan(f(:))) | any(isinf(f(:)))
    err = -2; 
  else
    err = checkvec(g, n);
    if  ~err,  err = checkvec(h, n); end
  end
end

function  err = checkvec(v,n)
sv = size(v); 
if  (min(sv) ~= 1) | ~isreal(v) | any(isnan(v(:))) | isinf(norm(v(:))) 
  err = -3; 
elseif  max(sv) ~= n,  err = -4;
else,                  err = 0; end