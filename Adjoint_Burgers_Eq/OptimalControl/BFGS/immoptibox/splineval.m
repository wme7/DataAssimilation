function  f = splineval(S, t)
%SPLINEVAL  Values of cubic spline.
% Call:  f = splineval(S, t)
%        f = splineval(t, S)
%
% Input
% S :  Struct defining the spline  s, see SPLINEFIT.
% t :  Vector with arguments for the spline.
%
% Output
% f :  Vector of the same type as  t , with  f(i) = s(t(i)) .

% Version 04.10.02  hbn(a)imm.dtu.dk

%  Check arguments
if  isstruct(t)  % reverse order
  temp = S;  S = t;  t = temp;
end
err = checksplxt(S,t,2,'SPLINEVAL');  % error stop if errors

% Sort arguments to compute from left to right
[t,js] = sort(t);   f = zeros(size(t));
pp = S.pp;   p = size(pp,2);   q = length(t);

%  Points to the left of the knot range
i = find(t <= pp(1,1));   usd = length(i);
if  usd
  v = t(i) - pp(1,1);
  f(js(i)) = pp(2,1) + v.*(pp(3,1) + pp(4,1)*v);
end
for  j = 1 : p-1
  i = find(t(usd+1:q) <= pp(1,j+1));   li = length(i);
  if  li    % Unused arguments
    v = t(usd+i) - pp(1,j);
    f(js(usd+i)) = pp(2,j) + v.*(pp(3,j) + v.*(pp(4,j) + pp(5,j)*v));
    usd = usd + li;
  end
end

if  usd < q  %  Points to the right of the knot range
  v = t(usd+1:q) - pp(1,p);
  f(js(usd+1:q)) = pp(2,p) + v.*(pp(3,p) + pp(4,p)*v);
end