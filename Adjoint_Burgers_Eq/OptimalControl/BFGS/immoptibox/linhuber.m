function  [X, info, perf] = linhuber(A,b, c, L, par, init)
%LINHUBER  Minimizer of the piecewise quadratic
%     f(x) = sum(phi(ri(x))) + c'*x + 0.5mu*||L*x||2^2 ,
% where ri(x) is the i'th element in r(x) = b - A*x and phi is the
% (possibly one-sided) Huber function with threshold gamma.
%  
% Call:    [X, info] = linhuber(A,b, c, L, par)
%          [X, info] = linhuber(A,b, c, L, par, init)
%          [X, info, perf] = linhuber(......)
%
% Input parameters
% A,b :  m*n matrix and m-vector, respectively.
% c   :  n-vector or empty, in which case the term c'*x is omitted.
% L   :  Matrix with n columns or empty.  In the case mu = par(2) > 0
%        an empty L is treated as L = I, the identity matrix.
% par :  Vector with one, two or three elements.
%        par(1) = gamma, Huber threshold.
%        par(2) = mu.  Default: mu = 0.
%        par(3) : Choice of Huber function,
%                 1 : one-sided,  rho(r) = 0 for r > 0,
%                 2 : one-sided,  all ri <= gamma are active
%                 Otherwise: standard Huber (default), ie 
%                            equations with  abs(ri)<=gamma  are active.
% init:  Choice of starting point,
%        init is a vector:  given x0,
%        init is a struct like info below: given active set and
%                                          factorization,
%        If nargin < 6, x0 is the least squares solution to Ax "=" b.   
%
% Output parameters
% X    :  If  perf  is present, then array, holding the iterates
%         columnwise.  Otherwise, computed minimizer of f.
% info :  Struct with information on the performance and computed solution.
%         Fields
%         pf :  Vector  [f(x)  ||g(x)||inf  #iterations #factorizations] ,
%               where x is the computed minimizer.
%               pf(1) = -Inf  indicates an unbounded problem.
%          S :  Struct with the Huber active set at the solution.  
%               Fields
%               s :  Huber sign vector,
%               A :  active set (indices of small residuals),
%               N :  inactive set (indices of large residuals),
%               L :  vector  [length(S.A), length(S.N)].
%          R :  n*n array with Cholesky factor of active set.
%          p :  Permutation vector used in the factorization.
% perf :  Iteration history.  Struct with fields
%          f :  values of  f(x) ,
%         ng :  values of  || g(x) ||inf
%         nA :  number of active equations.
%
% Method
% Essentially (without downdating of the factorization) as described 
% in H.B. Nielsen, "Computing a minimizer of a piecewise quadratic -
% Implementation", IMM-REP-1998-14, IMM, DTU.

% Version 10.10.05.  hbn(a)imm.dtu.dk

% Initialize
X = [];  S = [];  R = [];  p = [];
pf = [-Inf NaN 0 0];  perf = [];
[aux Htype] = lhaux(A,b,c,L,par);
[m n] = size(A);   kmax = 2*m;   fact = 0;
kcomp = 1;

if  nargin > 5  % given info on starting point
  if  isstruct(init)  % given active set and factorization
    kcomp = 2;
    S = init.S;
    [f g] = lhobj(zeros(n,1),b,S,c,L,A,aux);
    if  aux(3) == 0  % no L-contribution.  Reuse factorization
      [x u cf R p fact] = lhfacsolv(S,A,L,g,aux,S,init.R,init.p);  
    else
      [x u cf R p fact] = lhfacsolv(S,A,L,g,aux);  
    end
  else,  x = init(:); end  % given x0
else  % Start with least sqares solution to Ax "=" b
  x = pinv(A)*b;  fact = 1;
end

Trace = nargout > 2;
if  Trace
  X = zeros(n,kmax);  perf = zeros(3,kmax);
end 

% Iterate
k = 0;   stop = 0;

while  ~stop & k < kmax   
  Sp = S;
  k = k + 1;
  [r S tau] = lhsign(A,b,x,aux,Htype);
  [f g] = lhobj(x,r,S,c,L,A,aux);  ng = norm(g,inf);
  if  Trace,  X(:,k) = x;  perf(:,k) = [f; ng; S.L(1)]; end
  
  if  ng == 0 | ( k > kcomp & isequal(Sp,S) )
    stop = 1;  t = 0;
  else
    [h u lc R p rf] = lhfacsolv(S,A,L,g,aux,Sp,R,p);  fact = fact + rf;
    t = lhlines(S,r,u,lc,aux,Htype);
    if  t > 0,  x = x + t*h; 
    else,       stop = 1;    end
  end
end

if  t*ng == 0,  pf(1) = f; end
pf(2:4) = [ng k-1 fact];
info = struct('pf',pf, 'S',S, 'R',R, 'p',p);
if  Trace
  kk = 1 : k;  X = X(:,kk);  
  perf = struct('f',perf(1,kk), 'ng',perf(2,kk), 'nA',perf(3,kk)); 
else,       X = x; end

%%%%%%%%%%%%%%%%%%%%  Auxiliary functions  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [aux, Htype] = lhaux(A,b,c,L,par)
% Simple checks and set thresholds.  04.08.24
aux = [];
e = [A(:); b(:); c(:); L(:); par(:)];
v = isreal(e) & all(~isinf(e)) & all(~isnan(e));
if  ~v
  error('Indata contains non-real elements')
end
[m n] = size(A);  sb = size(b);
v = n > 0 & (min(sb) == 1 & max(sb) == m);
if  ~v
  disp('A must be a matrix and ')
  error(sprintf('b must be a vector of length %g',m))
end
sc = size(c);  sL = size(L);  tn = sprintf(' %g',n);
v = max(sc) == 0 | (min(sc) == 1 & max(sc) == n);
if  ~v
  error(['c must be empty or a vector of length' tn])
end
v = max(sL) == 0 | sL(2) == n;
if  ~v
  error(['L must be empty or a matrix with' tn ' columns'])
end
sp = size(par);
if  min(sp) ~= 1
  error('par must be a scalar or a vector')
end
if  max(sp) < 2,  par = [par 0]; end
if  ~all(par(1:2) >= 0)
  error('par(1:2) = [gamma  mu]  cannot be negative')
end

gamma = par(1);   mu = par(2);   gm = gamma*mu;
% Rounding error threshold
tau = 100*eps*norm(b,inf);   
if  par(1) < tau
  error(sprintf('gamma must be at least %9.2e',tau))
end

% Diverse norms and thresholds
nAi = norm(A,inf);  nA1 = norm(A,1);
sigma0 = 10*sqrt(eps * nA1*nAi);  sigma = sigma0;
if  gm > 0 
  sgm = sqrt(gm);
  if  isempty(L)  % L = I
    nL1 = 1;  nLi = 1;
    if  sgm >= sigma0, sigma = sgm; end
  else
    nL1 = norm(L,1);   nLi = norm(L,inf);
    if  gm*nL1*nLi >= sigma0^2
      nBB = (nA1 + sgm*nL1) * max(nAi, sgm*nLi);
      sigma = 10*sqrt(eps * nBB);
      sigma0 = sigma;
    end
  end
else
  nL1 = 0;  nLi = 0;
end  % mu-contribution

% aux = [gamma  tau  gamma*mu  sigma  sigma0   ...
%        ||A||inf  gamma*(||A||1 + mu*||L||1*||L||inf)  gamma*||c||inf]
aux = [gamma  tau  gm  sigma  sigma0 ...
    nAi  gamma*nA1 + gm*nL1*nLi  gamma*norm(c,inf)];
if  length(par) < 3 | abs(par(3) - 1.5) ~= 0.5,  Htype = 3;
else,                                            Htype = par(3); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [h,u,lc,R,p,refac] = lhfacsolv(S,A,L,g,aux,Sp,R,p)
% Factorize Hessian and find Newton step.  04.08.25
n = length(g);  
if  nargin < 6 | isempty(Sp),  refac = 1;
else  % possible update from previos factorization
  w = zeros(size(A,1),1);  w(S.A) = 1;  w(Sp.A) = w(Sp.A) - 1;
  if  any(w == -1),  refac = 1;  % downdates
  else
    refac = 0;   E = find(w);  % augmented active set
    if  ~isempty(E)  % update factorization
      [Q R p1] = qr([R; A(E,p)],0);  p = p(p1);
    end
  end
end

if  refac  % refactorize  
  d = aux(4)*ones(size(g));
  if  aux(3) > 0 & ~isempty(L)  % mu-contribution with L ~= I
    C = [A(S.A,:); sqrt(aux(3))*L; diag(d)];
  else
    C = [A(S.A,:); diag(d)];
  end
  [Q R p] = qr(C,0);
end

% Check for almost rank deficient
rd = any(abs(diag(R)) < 10*aux(5));
% Get solution to  Hx = -g
g = -g(:);   x(p) = R \ (g(p)' / R)';  x = x(:);
if  rd  % Check consistency
  v(p) = R \ (aux(4)^2*x(p)' / R)';   v = v(:);  h = x + v;
  if  S.L(1) > 0,  gR = A(S.A,:)' * (A(S.A,:)*h);
  else,            gR = zeros(n,1);               end
  if  aux(3) > 0  % mu-contribution
    if  isempty(L),  gR = gR + aux(3)*h;
    else,            gR = gR + aux(3)*(L' * (L*h)); end
  end
  if  norm(g - gR) > sqrt(eps)*norm(g)  % inconsistent (?)
    h = v / norm(v);
  end
else  % full rank.  iterative refinement step
  if  S.L(1) > 0,  Hx = A(S.A,:)' * (A(S.A,:)*x);
  else,            Hx = zeros(n,1);               end
  if  aux(3) > 0  % mu-contribution
    if  isempty(L),  Hx = Hx + aux(3)*x;
    else,            Hx = Hx + aux(3)*(L' * (L*x)); end
  end
  v(p) = R \ ((g(p) - Hx(p))' / R)';   h = x + v(:);
end

% u = A*h, modified for 'zero' elements
u = A*h;
z = find(abs(u) <= eps * aux(4) * norm(h,inf) );
u(z) = 0;

% Offset values for line search coefficients
if  aux(3) > 0  % mu-contribution
  if  isempty(L),  c2 = aux(3)*norm(h)^2;
  else,            c2 = aux(3)*norm(L*h)^2; end
else
  c2 = 0;
end
lc = [dot(g,h)  c2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [kv, lc] = lhkinks(S,r,u,lc,aux,Htype)
% Get kink values and finish initialization of line search coeffs.  04.08.23
gam = aux(1);  tau = aux(2);  
if  nargin < 6 | abs(Htype - 1.5) ~= 0.5,  Htype = 3; end
thr = [-1 0; -inf 1; -1 1];  thr = gam*thr(Htype,:);  % Huber thresholds
kv = zeros(2*length(r),2);  k = 0;
if  ~isempty(S.A)
  % Remove leaving 'maybes'
  act = S.A;   w = ones(size(act));
  o = find((u(act) > 0 & abs(r(act) - thr(1)) <= tau) | ...
    (u(act) < 0 & abs(r(act) - thr(2)) <= tau));
  w(o) = 0;   act = act(find(w));
  if  ~isempty(act)  % remaining actives
    lc(2) = lc(2) + norm(u(act))^2;  % update initial line search coeffs.
    n = find(u(act) < 0);   p = find(u(act) > 0);  
    ii = act([p; n])';   li = length(ii);
    if  li  % get leave values
      d = [repmat(thr(1),length(p),1); repmat(thr(2),length(n),1)];
      kv(k+[1:li],:) = [(r(ii) - d)./u(ii) -ii];   k = k + li;
    end
  end
end  % actives

% Enter-leave kinks
i = find(u(S.N) .* r(S.N) > 0);   E = S.N(i);
if  Htype == 1  % add neglected components
  E = [E(:); find(r > tau & u > 0)]';
end
if  ~isempty(E)
  i = E(find(u(E) > 0));  j = E(find(u(E) < 0));  
  ii = [i  j]';   li = length(ii);
  d = [repmat(thr([2 1]),length(i),1); repmat(thr,length(j),1)];
  kv(k+[1:li],:) = [(r(ii) - d(:,1))./u(ii) ii];   k = k + li;
  kv(k+[1:li],:) = [(r(ii) - d(:,2))./u(ii) -ii];   k = k + li;
end

% sort by kink values
[sk s] = sort(kv(1:k,1));
% if  k > 2 & kv(s(k),2) < 0,  s = s(1:k-1); end % remove last leave
kv = kv(s,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  t = lhlines(S,r,u,lc,aux,Htype)
% Piecewise linear line search.  04.08.23
[kv lc] = lhkinks(S,r,u,lc,aux,Htype); 
c1 = lc(1);  c2 = lc(2);
tj = 0;  j = 0;  nk = size(kv,1);  fd = 0;
while  ~fd & j < nk
  j = j+1;  t1 = tj;  tj = kv(j,1);
  dt = tj - t1;  cj = c2*dt - c1;
  if  c2 & (cj >= 0)  % passed minimum
    t = t1 + c1/c2;   fd = 1;
  else  % update coefficients
    c1 = -cj;  ii = kv(j,2);
    c2 = c2 + sign(ii)*u(abs(ii))^2;
  end
end
if  fd,           return
elseif  nk > 1,   t = mean(kv(nk-1:nk,1));
elseif  nk == 1,  t = kv(1,1) + aux(2)/abs(u(abs(kv(1,2))));
elseif  c2,       t = c1/c2;
else,             t = -1;  end  % error return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [f,g] = lhobj(x,r,S,c,L,A,aux)
% Value and gradient of objective function.  04.08.24

% Compute function
fA = r(S.A);  fN = r(S.N);  gamma = aux(1);
f = norm(fA)^2/(2*gamma) + norm(fN,1) - .5*gamma*S.L(2);
if  ~isempty(c),  f = f + dot(c,x); end
if  aux(3) > 0
  if  isempty(L),  Lx = x; else,  Lx = L*x; end
  f = f + aux(3)/(2*gamma)*norm(Lx)^2;
end

if  nargout > 1  % gamma*gradient
  fv = -gamma*S.s';  fv(S.A) = -fA;  g = A'*fv;
  if  ~isempty(c),  g = g + gamma*c; end
  if  aux(3) > 0
    if  isempty(L),  g = g + aux(3)*x;
    else,            g = g + L'*(L*(aux(3)*x)); end
  end
  % Threshold for rounding errors
  thr = 100*eps*(aux(7)*norm(fv,inf) + aux(7)*norm(x,inf) + aux(8));
  if  norm(g,inf) <= thr,  g = zeros(size(g)); end
end  % gradient

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [r,S,tau] = lhsign(A,b,x,aux,Htype)
% Residual and Huber sign vector etc

% Rounding error and effective thresholds
tau = max(aux(2), 10*eps*aux(6)*norm(x,inf));
thr = aux(1) + tau;

% Residual and sign vector
r = b - A*x;   s = sign(r);

if  nargin < 5 | abs(Htype - 1.5) ~= 0.5  % simple Huber function
  w = find(abs(r) <= thr);    s(w) = 0;   
elseif  Htype == 1  % Neglect positive contributions
  w = find(-thr <= r & r <= tau);    s(w) = 0;
  p = find(s > 0);   s(p) = 0;
else                % Negative contributions are active
  w = find(r <= thr);   s(w) = 0;
end
N = find(s);

% Return result
S = struct('s',s', 'A',w', 'N',N', 'L',[length(w)  length(N)]);