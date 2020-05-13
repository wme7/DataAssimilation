function  S = splinefit(tyw, x, D)
%SPLINEFIT  Weighted least squares fit with cubic spline  
% Call    S = splinefit(tyw, x)
%         S = splinefit(tyw, x, D)
%
% Input
% tyw :  Data points and weights.  Array with 2 or 3 columns,
%        tyw(:,1) :  Abscissas.
%        tyw(:,2) :  Ordinatates.
%        tyw(:,3) :  Weights.  If  tyw  holds less than 3 columns, then
%                    all weights are set to 1.
% x   :  Knots.  Must satisfy  x(1) < x(2) <= ... <= x(end-1) < x(end) .
% D   :  If present and nonempty, it must either be a 
%        scalar:  Periodic spline if D=1, otherwise a natural spline,
%        or a 2*4 matrix: The fit is computed with respect to
%        boundary conditions
%            D(k,1:3)*[s(zk); s'(zk); s"(zk)] = D(k,4)  for  k=1,2
%        with  z1=x(1), z2 = x(end).
%
% Output
% S :  Struct representing the spline.  Fields
%      fail :  Performance indicator.  The other fields are only
%              significant if S.fail = 0.
%              fail = 0 :  No problems.
%                     1 :  x is not a real valued vector of 
%                          length at least 2.
%                     2 :  Knots are not in increasing order.
%                     3 :  Some data abscissa is not real or is
%                          outside [x(1), x(end)].
%                     4 :  No data ordinates are given.
%                     5 :  Too few active data points (points with
%                          strictly positive weights).
%                     6 :  Schoenberg-Whitney condition is not satisfied.
%                     7 :  When given, D should be either a 
%                          scalar or a 2*4 matrix.
%      x    :  Knots.
%      c    :  Coefficients in B-spline representation of  s .
%      pp   :  Piecewise polynomial representation of  s .
%      sdy  :  Estimated standard deviation of data.
%      sdc  :  Estimated standard deviation of c.

% Version 10.10.05.  hbn(a)imm.dtu.dk

% Initialize
S = struct('fail',0, 'x',x, 'c',NaN, 'pp',NaN, 'sdy',NaN, 'sdc',NaN);

% Check knots and data 
[S.fail, tyw, ma] = checksplxt(x,tyw,1, 'SPLINEFIT');
if  S.fail, return, end

% Extended knot set 
x = x(:)';   n = length(x) - 1;  % knots as row vector.  n intervals
xx = [x(1)*[1 1] x x(end)*[1 1]];

%  Possible boundary conditions
period = 0;  ns = n+3;
con = zeros(1,2);    alfa = zeros(1,6);
if  nargin > 2 & ~isempty(D)
  sd = size(D);
  if  all(sd == 1)             % scalar input
    if  D == 1                 % periodic spline
      period = 1;  ns = n;
    else                       % natural spline
      D = [0 0 1 0; 0 0 1 0];  sd = size(D);
    end
  end
  if  ~period
    if  any(sd - [2 4]),  S.fail = 7;  return,  end
    % general boundary conditions
    B = splbsd(x(1), xx(1:6));
    alfal = D(1,1:3)*B(:,1:3);
    if  alfal(1)               % Active constraint
      alfa(1:3) = [-alfal(2:3) D(1,4)]/alfal(1);   con(1) = 1;  
    end
    B = splbsd(x(n+1), xx(n:n+5));
    alfar = D(2,1:3)*B(:,2:4);
    if  alfar(3)               % Active constraint
      alfa(4:6) = [-alfar(1:2) D(2,4)]/alfar(3);   con(2) = 1;
    end
    ns = ns - sum(con);
  end
end  %  Boundary conditions

% Check number of active points
if  ma < ns,  S.fail = 5;  return, end

% Get elements of matrix and right hand side
N = zeros(ma,5);  j1 = ones(ma,1);
i2 = 0;  % currently last point
for  j = 1 : n
  i = find(tyw(i2+1:ma,1) <= x(j+1));   li = length(i);
  if  li    %  Contributions
    i1 = i2+1;  i2 = i2+li;  ii = i1:i2;
    N(ii,:) = repmat(tyw(ii,3),1,5) .* ...
      [splbsv(tyw(ii,1),xx(j:j+5)) tyw(ii,2)];
    j1(ii) = j;
  end
end

% Store in sparse matrix and rhs
nz = 4*ma;
ii = reshape(repmat(1:ma,4,1), nz,1);
j = [j1 j1+1 j1+2 j1+3];   jj = reshape(j',nz,1);
F = sparse(ii,jj,reshape(N(:,1:4)',nz,1),ma,n+3,nz);   b = N(:,5);

if  period  % Reduce for periodicity conditions  
  B1 = splbsd(x(1), xx(1:6));   Bm = splbsd(x(end),xx(end-5:end));
  B = Bm(:,2:4)\B1(:,1:3);
  F(:,1:3) = F(:,1:3) + F(:,n+[1:3])*B;  F = F(:,1:n);
end
if  con(2)  % reduce at right hand end
  ii = find(F(:,n+3));
  if  ~isempty(ii),  b(ii) = b(ii) - alfa(6)*F(ii,n+3); end
  F = [F(:,1:n) F(:,n+1:n+2)+F(:,n+3)*alfa(4:5)];
end
if  con(1)  % reduce at left hand end
  ii = find(F(:,1));
  if  ~isempty(ii),  b(ii) = b(ii) - alfa(3)*F(ii,1); end
  F = [F(:,2:3)+F(:,1)*alfa(1:2) F(:,4:end)];
end

% Get solution with check of full rank
scol = sum(F.^2);
if  any(scol == 0)
  S.fail = 6;  % Schoenberg-Whitney not satisfied
  return
end
R = qr([F b],0);
reldR = diag(R(1:ns,1:ns)).' ./ scol;
if  any(abs(reldR) < sqrt(eps))
  S.fail = 6;  % Schoenberg-Whitney not satisfied
  return
end

% Data variance
if  ma > ns,   vy = R(end,end)^2/(ma - ns);
else,  vy = 0; end

% Get solution and coefficient variance
c = ones(n+3,1);  vc = vy*c;
if      period,  ic = ns;
elseif  con(2),  ic = n+2;  
else,   ic = n+3;  end
dj = ic - ns;
j = 0;  % number of upper triangular elements in row
for  i = ns : -1 : 1
  jj = i + (1:j);  
  c(ic) = (R(i,end) - R(i,jj)*c(jj+dj))/R(i,i);
  vc(ic) = (vy + R(i,jj).^2*vc(jj+dj))/R(i,i)^2;
  ic = ic -1;  j = min(j+1,3);
end
if  ns < n+3  % get remaining coefficients
  if  period
    c(n+[1:3]) = B*c(1:3);
    vc(n+[1:3]) = B.^2 * vc(1:3);
  elseif  con(2)
    c(n+3) = alfa(4:6)*c(n+(1:3));
    vc(n+3) = alfa(4:6).^2 * vc(n+(1:3));
  end
  if  con(1)
    jj = [2 3 1];
    c(1) = alfa(1:3)*c(jj);
    vc(1) = alfa(1:3).^2 * vc(jj);
  end
end

% Return results
S.x = x;   S.c = c.';   
S.pp = splmpp(x,c);
if  period,  S.pp(2:end,end) = S.pp(2:end,1); end
S.sdy = sqrt(vy);   S.sdc = sqrt(vc)';

% ==========  auxiliary functions  =====================================

function  N = splbsd(t,z)
% Values and first two derivatives of the four nonzero B-splines
% at (scalar)  t,   z3 <= t <= z4
M = zeros(3,5);   N = zeros(3,4);
M(1,2) = 1/(z(4) - z(3));
for  r = 1 : 2
  s = 0 : r;   jj = 3 - r;
  M(r+1,2+s) = ((t - z(jj+s)).*M(r,1+s) + ...
    (z(4+s) - t).*M(r,2+s)) ./ (z(4+s) - z(jj+s));
end
N(1,1:3) = (z(4:6) - t) .* M(3,2:4);
N(1,2:4) = N(1,2:4) + (t - z(1:3)) .* M(3,2:4);
%  Differentiate
N(2,:) = -3*diff(M(3,:));
M(2,2:4) = diff(M(2,1:4)) ./ (z(4:6) - z(1:3));
N(3,:) = 6*diff(M(2,:));

% ======================================================================

function  N = splbsv(t,z)
% Values of the four nonzero B-splines at (vector)  t,  
% z3 <= t <= z4
q = length(t);   N = zeros(q,4);
z3 = z(3);   z4 = z(4);   hj = z4 - z3;
u = t - z3;   v = z4 - t;
N(:,1) = v/(hj*(z4 - z(2)));
N(:,2) = u/(hj*(z(5) - z3));    % M1
N(:,3) = u.*N(:,2)/(z(6) - z3);
N(:,2) = ((t - z(2)).*N(:,1) + (z(5)-t).*N(:,2))/(z(5)-z(2));
N(:,1) = v.*N(:,1)/(z(4)-z(1));    % M2
N(:,4) = u.*N(:,3);
N(:,3) = (t - z(2)).*N(:,2) + (z(6) - t).*N(:,3);
N(:,2) = (t - z(1)).*N(:,1) + (z(5) - t).*N(:,2);
N(:,1) = v.*N(:,1);

% ======================================================================

function  pp = splmpp(x,c)
% Given knots  x(1:n+1) and B-spline coefficients c(1:n+3).
% Compute piecewise polynomial representation.
n = length(x) - 1;   pp = zeros(5,n+1);   p = 0;
xx = [x(1) x(1) x(:)' x(n+1) x(n+1)];
for  j = 1 : n
  z = xx(j:j+5);
  if  x(j+1) > x(j)  % Non empty interval
    p = p+1;   pp(1,p) = x(j);
    if  z(2) == x(j)  % After empty interval
      B = splbsd(x(j),z);   q = 1;
    else,  q = 2;  end
    pp(2:4,p) = B(:,q:q+2)*c(j:j+2);   pp(4,p) = .5*pp(4,p);
    B = splbsd(x(j+1),z); 
    d2s = .5*(B(3,2:4)*c(j+1:j+3));
    pp(5,p) = (d2s - pp(4,p))/(3*(x(j+1) - x(j)));
  end  % Non empty interval
end  % j
p = p+1;   pp(1,p) = x(n+1);    pp(2,p) = c(n+3);
pp(3,p) = B(2,3:4)*c(n+2:n+3);  pp(4,p) = d2s;
pp = pp(:,1:p); 