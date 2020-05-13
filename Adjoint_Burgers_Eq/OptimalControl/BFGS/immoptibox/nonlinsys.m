function  [X, info, perf] = nonlinsys(fun, x0, opts, B0, varargin)
%NONLINSYS Powell's dog-leg method for non-linear system of equations.
% Find an n-vector xm, such that f_i(xm) = 0 , i=1,...,n.
% The vector function  f(x) (with components f_i, i=1,...,n) must be 
% given by a MATLAB function with declaration
%            function  f = fun(x, p1,p2,...)
% % % p1,p2,... are parameters of the function.  
%  
% Call
%    [X, info] = nonlinsys(fun, x0)
%    [X, info] = nonlinsys(fun, x0, opts)
%    [X, info] = nonlinsys(fun, x0, opts, B0, p1,p2,...)
%    [X, info, perf] = nonlinsys(.....)
%
% Input parameters
% fun  :  Handle to the function.
% x0   :  Starting guess for  xm .
% opts :  Vector with five elements:
%         opts(1)    initial trust region radius Delta
%         opts(2:4)  used in stopping criteria: 
%             ||f||inf <= opts(2)                      or
%             ||dx||2 <= opts(3)*(opts(3) + ||x||2)    or
%             no. of iteration steps exceeds  opts(4) . 
%         opts(5)  "relative" step length for difference approximations.
%         Default  opts = [0.1(1+||x0||) 1e-6 1e-8 100 1e-6]
%         If the input opts has less than 5 elements, it is
%         augmented by the default values. 
% B0   :  Initial approximation to the Jacobian of f.
%         If  B0 is not given, a forward difference approximation
%         to it is used.
% p1,p2,..  are passed dirctly to the function FUN .    
%
% Output parameters
% X    :  If  perf  is present, then array, holding the iterates
%         columnwise.  Otherwise, computed solution vector.
% info :  Performance information, vector with 6 elements:
%         info(1:3) = final values of 
%             [||f(x)||inf  ||dx||2  Delta] 
%         info(4:5) = no. of iteration steps and function avaluations
%         info(6) = 1 :  Stopped by small f-vector
%                   2 :  Stopped by small x-step
%                   3 :  No. of iteration steps exceeded
%                  -1 :  x is not a real valued vector
%                  -2 :  f is not a real valued column vector
%                  -3 :  Dimension mismatch in x, f, B0
%                  -4 :  Maybe started at a saddle point
%                  -5 :  Overflow during computation 
% perf :  Array, holding 
%         perf(1,:) = values of  || f(x) ||inf
%         perf(2,:) = Delta-values.

% Version 04.04.10.  hbn(a)imm.dtu.dk

% Check parameters and function call
if  nargin < 2
  [stop,x,f,opts,B,D,info,n] = checkinput(fun, [], [], []);
elseif  nargin < 3
  [stop,x,f,opts,B,D,info,n] = checkinput(fun, x0, [], []);
elseif  nargin < 4
  [stop,x,f,opts,B,D,info,n] = checkinput(fun, x0, opts, []);
else
  [stop,x,f,opts,B,D,info,n] = checkinput(fun, x0, opts, B0, varargin{:});
end
if  stop
  perf = [];  X = x;  return
end

% Finish initialization
Delta = opts(1);    kmax = opts(4);
Trace = nargout > 2;
if  Trace
  X = repmat(x,1,kmax+1);
  perf = repmat([norm(f,inf); Delta],1,kmax+1);
end 
k = 1;   nx = norm(x);  nf = norm(f,inf);  F = .5*norm(f)^2;
ku = 0;  refac = 0;  % For "extra" updates
newD = 0;  % Signifies whether D should be recomputed
nstep = 0;
% Iterate
while   ~stop 
  if      isinf(nf),                        stop = -5;
  elseif  nf <= opts(2),                    stop = 1;  
  elseif  Delta <= opts(3)*(opts(3) + nx),  stop = 2; 
  elseif  k >= kmax,                        stop = 3;
  else
    % Newton step
    if  newD,  [stop D hN] = getDandh(B,f);
    else,      hN = D*(-f);  end
    nhN = norm(hN);
    if  nhN <= opts(3)*(opts(3) + nx),  stop = 2;
    elseif  ~stop
      if  nhN > Delta  % Include gradient
        g = B'*(-f);  ng = norm(g);  alpha = (ng / norm(B*g))^2;
        gn = alpha*g;  ngn = alpha*ng;
        if  ngn >= Delta
          h = (Delta/ngn) * gn;
        else  %  Dog leg
          b = hN - gn;  bb = b'*b;  gb = gn'*b;
          c = (Delta + ngn)*(Delta - ngn);
          if  gb > 0
            beta = c / (gb + sqrt(gb^2 + c * bb));
          else
            beta = (sqrt(gb^2 + c * bb) - gb)/bb;
          end
          h = gn + beta*b;
        end
      else
        h = hN;  
      end
    end
  end
  if  ~stop
    ku = mod(ku,n) + 1;  
    if  abs(h(ku)) < .8*norm(h)  % extra step for updating B and D
      xu = x;
      if  x(ku) == 0,  xu(ku) = opts(5)^2;
      else,            xu(ku) = x(ku) + opts(5)*abs(x(ku)); end
      [stop,B,D,Fu,fu,hu,nhu,newD] = updBD(fun,x,xu,f,B,D,0,varargin{:});
      info(5) = info(5)+1;
    end
    if  ~stop
      xnew = x + h;      
      [stop,B,D,Fn,fn,h,nstep,newD] = updBD(fun,x,xnew,f,B,D,newD,varargin{:});
      info(5) = info(5)+1;
      if  ~stop
        % Check trust region
        dF = F - Fn;
        if  dF > 0                        % Update x  
          dL = F - .5*norm(f + B*h)^2; 
          x = xnew;   F = Fn;  f = fn;  
          nf = norm(f,inf);   nx = norm(x);
          if  dL > 0,  rho = dF / dL;  else,  rho = 0; end
        else,  rho = 0;   end
        if      rho > .75,  Delta = max(Delta, 3*nstep);
        elseif  rho < .25,  Delta = Delta/2; end
        k = k + 1;
        if  Trace,  X(:,k) = x;   perf(:,k) = [nf Delta]'; end
      end
    end
  end   
end
%  Set return values
if  Trace
  X = X(:,1:k);   perf = perf(:,1:k);
else,  X = x;  end
if  stop < 0,  nf = NaN; end
info = [nf  nstep  Delta  k-1  info(5) stop];

%%%%%%%%%%%%%%%%%%%%  Auxiliary functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [stop,x,f,opts,B,D,info,n] = checkinput(fun, x0, opts, B0, varargin)
% Check indata
info = [NaN zeros(1,5)];  f = NaN;  B = NaN;  D = NaN;  n = NaN;
stop = 0;
if  isempty(x0),  stop = -1;  x = [];
else
  [stop x n] = checkx(x0);   
  if  ~stop
    [stop F f] = checkfJ(fun,x0,varargin{:});  
    info([1 5]) = [norm(f,inf) 1];
    nf = norm(f, inf);
    if  ~stop & length(f) ~= n,  stop = -3; end
    if  ~stop
      %  Finish initialization
      opts = checkopts(opts, [.1*(1+norm(x,inf)) 1e-6 1e-8 100 1e-6]); 
      % Jacobian
      sB = size(B0); 
      if  sum(sB) == 0  % placeholder
        [stop B] = Dapprox(fun,x,opts(5),f,varargin{:});  
        info(5) = info(5) + n;
      elseif  any(sB ~= n),  stop = -3;
      else,                  B = B0;   end
      % Check gradient 
      if  ~stop
        g = B'*(-f);   ng = norm(g,inf);  
        if  isinf(ng),  stop = -5; end 
      end
      % Get initial inverse Jacobian and check D*f
      if  ~stop
        [stop D hN] = getDandh(B,f);
      end
    end
  end
end
if  stop,  info(6) = stop; end

function  [stop, D, hN] = getDandh(B,f)
% Get inverse of approximate Jacobian, and check for overflow
[U S V] = svd(B);
s = diag(S);  i = find(s > 100*eps*s(1));
if  isempty(i),  stop = -4;  D = NaN;  hN = NaN;
else
  D = V(:,i) * diag(1./s(i)) * U(:,i)';
  hN = D*(-f);
  if  isinf(norm(hN)),  stop = -5; else,  stop = 0; end
end

function  [stop,B,D,Fn,fn,h,nh,newD] = updBD(fun,x,xn,f,B,D,newD,varargin)
% Evaluate at new point and update B and D
[stop Fn fn] = checkfJ(fun,xn,varargin{:});  
if  ~stop
  h = xn - x;  nh = norm(h);  y = fn - f;
  B = B + ((y - B*h)/nh^2) * h';
  if  ~newD
    Dy = D*y;   hDy = dot(h,Dy);
    if  abs(hDy) < sqrt(eps)*nh, newD = 1;
    else
      D = D + ((h - Dy)/hDy)*(h'*D);
    end
  end
end