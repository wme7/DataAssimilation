function  [Z, c, info, perf] = mexpfit(tyw, z0, opts)
%MEXPFIT  Least squares fit with a sum of decaying exponentials:  Find  
% zm = argmin{f(z)} , where  z  is a p-vector with positive elements, 
%     f(z) = ½||W*r(z)||^2 ,   r(z) = y - F(z)*c(z) ,
% where  F(i,:) = [exp(-z1*ti)  ...  exp(-zp*ti)]      or
%        F(i,:) = [exp(-z1*ti)  ...  exp(-zp*ti)  1] ,
% depending on the option 'const', and  c(z) is the least squares 
% solution to  W*F*c "=" W*y .
%
% Call    [Z, c, info] = mexpfit(tyw, z0)
%         [Z, c, info] = mexpfit(tyw, z0, opts)
%         [Z, c, info, perf] = mexpfit(...)
%
% Input parameters
% tyw :   Data points and weights.  Array with 2 or 3 columns,
%         tyw(:,1:2) :  Abscissas and ordinates of data points.
%         tyw(:,3)   :  Weights.  If  tyw  has less than 3 columns, then
%                       all weights are set to 1.
% z0   :  Starting guess for  z .
% opts :  Either a struct with fields 'const', 'tau', 'tolg', 'tolx', 
%         and 'maxeval', or a vector with the values of these options,
%         opts = [const  tau  tolg  tolx  maxeval].
%         const :  If positive then there is a constant term.
%         tau   :  Used in starting value for Marquardt parameter for the
%                  optimization,  mu = tau * max( [J(z0)'*J(z0)]_(ii) ) . 
%         The other options are used in stopping criteria:
%             ||g||_inf <= tolg                            or 
%             ||dz||_2  <= tolx*(tolx + ||z||_2)           or
%             no. of function evaluations exceeds  maxeval .
%         Default  const = 0,  tau = 1e-3,  tolg = 1e-6,  tolx = 1e-10,  
%         maxeval = 100.  If the input opts has less than 5 elements, it
%         is augmented by the default values.  Also, zeros and negative 
%         elements are replaced by the default values.
%
% Output parameters
% Z    :  If  perf  is present, then an array, holding the z-iterates
%         columnwise.  Otherwise, computed solution vector zm.
% c    :  If  info.fail = 0, then c holds coefficients [c1, ..., cq], where 
%         opts.const > 0 :  q = p+1;  cq is the constant term
%              otherwise :  q = p .
% info :  Struct with information on performance.  Fields
%         fail:  Successful run if fail = 0.  Otherwise info.msg tells
%                what went wrong.
%          msg:  Text message about performance.
%         vals:  Vector with final values of  f(z), ||g(z)||inf, ||dz||2
%          its:  Vector with number of iterations and function evaluations.
% perf :  Struct with fields
%          f :  values of  f(zk) ,
%         ng :  values of  || g(zk) ||_inf ,
%         mu :  values of Marquardt damping parameter.

% Version 10.10.29.  hbn(a)imm.dtu.dk

%  Initialise with check of input
Z = [];  c = [];  perf = [];
if  nargin < 3,  opts = []; end
[z info opts nc t y W] = init(tyw, z0, opts);
if  info.fail,  return, end
p = length(z);

if  p == 1 & nc == 1 % simple case
  [Z c f] = simple(t, y, W);
  if  ischar(c),  info.fail = 1;  info.msg = c;  
  else,           info.vals(1) = f;  info.its = [0 1];   end
  return
end

% General case.  Remaining options and initial values
tau = opts(2);    tolg = opts(3);  tolx = opts(4);  maxeval = opts(5); 
[r J c] = func(z, t,y,W, nc);  nv1 = 1; 
A = J'*J;   g = J'*r;   f = (r'*r)/2;   ng = norm(g,inf);
mu = tau * max(diag(A)); 
  
Trace = nargout > 3;
if  Trace
  Z = z*ones(1,maxeval);   perf = [f; ng; mu]*ones(1,maxeval);
end

nu = 2;  mok = 0;  nh = 0;  n = size(A,1);  stop = 0;

while   ~stop
  if  ng <= tolg,  stop = 1;
  else
    hh = (A + mu*eye(n))\(-g);
    nh = norm(hh);   nz = tolx + norm(z);
    if      nh <= tolx*nz,  stop = 2;
    elseif  nh >= nz/eps,   stop = 4; end    % Almost singular ?
  end
  if  ~stop
    % Control h to guarantee znew > 0
    if  p == 1
      hneg = -0.75*z;  hpos = z;
    else
      hneg = ([z(2:p); -z(p)] - z)/3;  
      hpos = ([2*z(1); z(1:p-1)] - z)/3; 
    end
    hmod = max( min(hh, hpos), hneg);
    znew = z + hmod;   h = znew - z;   
    dr = J*h;  dL = -(dr' * (2*r + dr));  % twice linear gain
    [rn Jn c] = func(znew, t,y,W, nc);  nv1 = nv1 + 1;
    if  ischar(c),  info.fail = 1;  info.msg = c;  return, end
    fn = (rn'*rn)/2;   df = (r - rn)' * (r + rn);  % twice actual gain
    if  (dL > 0) & (df > 0)
      mok = mok + 1;
      z = znew;   f = fn;  J = Jn;  r = rn;
      A = J'*J;   g = J'*r;   ng = norm(g,inf);
      if  norm(hh - hmod)  % step was reduced.  Increase mu
        mu = mu*nu;  nu = 2*nu;
      else
        mu = mu * max(1/3, 1 - (2*df/dL - 1)^3);   nu = 2;
      end
      if  Trace, Z(:,mok+1) = z;  perf(:,mok+1) = [f; ng; mu]; end
    else    % Marquardt fail
      mu = mu*nu;  nu = 2*nu;
      if  mok > n & nu > 8, stop = 5; end
    end
    if  nv1 >= maxeval,  stop = 3; end
  end
end
%  Set return values
if  Trace
  jj = 1 : (mok+1);  Z = Z(:,jj);   
  perf = struct('f',perf(1,jj), 'ng',perf(2,jj), 'mu',perf(3,jj));
else,  Z = z;  end
if  stop == 1,      info.msg = 'Iterations stopped by a small gradient';
elseif  stop == 2,  info.msg = 'Iterations stopped by a small z-step';
elseif  stop == 3,  info.msg = 'Iterations stopped by  maxeval';
elseif  stop == 4,  info.msg = 'Iterations stopped by extreme step';
elseif  stop == 5,  info.msg = 'Iterations stopped by stalling';       end
info.vals = [f ng nh];
info.its = [mok nv1]; 

% ==========  auxiliary functions  =================================

function  [z, info, opts, nc, t,y,W] = init(tyw, z0, opts)
Z = NaN;  t = [];  y = [];  W = []; 
info = struct('fail',0, 'msg','', 'vals',[NaN NaN NaN], 'its',[0 0 0]);
sz = size(z0);
if  ~isreal(z0) | min(sz) ~= 1
  info.fail = 1;  
  info.msg = 'z must be a vector with real elements';  return
end
if  any(z0 <= 0)
  info.fail = 1;  
  info.msg = 'z must must have strictly positive elements';  return
end
opts = checkopts('mexpfit', opts);
[m q] = size(tyw);  p = length(z0);
if  opts(1),  nc = p+1;  else,  nc = p; end
if  q > 2  % remove points with non-positive weight
  i = find(tyw(:,3) > 0);  ma = length(i);
  if  ma < m,  tyw = tyw(i,:);  m = ma; end
else,  tyw = [tyw ones(m,1)]; end
if  p + nc > m
  info.fail = 1;  info.msg = 'Underdetermined problem';  return
end
t = tyw(:,1);  y = tyw(:,2);  W = tyw(:,3)*ones(1,p);
z = sort(z0(:),1,'descend');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [r, J, c] = func(z, t,y,W, nc)
% Function for exponential fit

p = length(z);  m = length(t);
E = W .* exp(t * (-z'));  y = W(:,1).*y;
if  nc > p,  F = [E W(:,1)];  else,  F = E; end
[Q R] = qr(F, 0);  rc = rcond(R);  pp = size(R,2);
if  rc < max(pp,10) * eps
  r = [];  J = [];
  c = 'Rank deficient system';  return
 elseif  rc < sqrt(eps)  % use iterative refinement 
  md = max(abs(diag(R))) * sqrt(eps);
  RR = qr([R; md*eye(pp)]);  RR = triu(RR(1:pp,:));
  c = RR \ (y' * Q)';
  r = y - F*c;
  c = c + RR \ (r' * Q)';
else
  md = 0;
  c = R \ (y' * Q)';
end
r = y - F*c;  % residual.  Now get Jacobian - without iterative refinement
tE = (t * ones(1,p)) .* E;  rho = (r' * tE)';
tEc = tE .* (ones(m,1) * c(1:p)');
if  nc > p,  B = [diag(rho); zeros(1,p)];  else,  B = diag(rho); end
J = tEc - Q * (R' \ (F'*tEc - B));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [z, c, f] = simple(t, y, w)
% One exponential only.  Use logaritmic transformation

% Ignore negative ordinates
m = length(y);  a = find(y > 0);  ma = length(a);
if  ma < 2
  c = 'Simple case.  Less than 2 active points';
  f = NaN;  return
end

if  ma < m,  t = t(a);  y = y(a);  w = w(a); end
wt = w .* y;  % transformed weights
x = [wt  wt.*t] \ (wt .* log(y));
z = -x(2);  c = exp(x(1));
r = w .* (y - c * exp(-z*t)); 
f = 0.5 * norm(r)^2;