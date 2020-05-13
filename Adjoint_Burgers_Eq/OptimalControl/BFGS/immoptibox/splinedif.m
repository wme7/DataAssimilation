function  f = splinedif(S, t, d)
%SPLINEDIF  Derivative of cubic spline.
% Call:    f = splinedif(S, t)
%          f = splinedif(S, t, d)
%
% Input
% S :  Struct defining the spline  s, see SPLINEFIT.
% t :  Vector with arguments. 
% d :  Order of differantiation.  Default d = 1.
%
% Output
% f :  Vector of the same type as  t , with  f(i) = s^(d)(t(i)) .

% Version 04.04.23.  hbn(a)imm.dtu.dk

% Check arguments
err = checksplxt(S,t,3,'SPLINEDIF');  % error stop if errors

% Get differentiation order
if  nargin < 3 | max(size(d)) ~= 1,  d = 1;
elseif  ~isreal(d) | isnan(d) | isinf(d) | d ~= round(d) | d < 0
  d = 1; 
end

if  d > 3,  f = zeros(size(t));
else
  if  d > 0
    switch  d
      case 1,  A = [0 1 0 0; 0 0 2 0; 0 0 0 3; 0 0 0 0];
      case 2,  A = [0 0 2 0; 0 0 0 6; 0 0 0 0; 0 0 0 0];
      case 3,  A = [0 0 0 6; 0 0 0 0; 0 0 0 0; 0 0 0 0];
    end
    S.pp(2:5,:) = A*S.pp(2:5,:);
  end
  f = splineval(S,t);
end