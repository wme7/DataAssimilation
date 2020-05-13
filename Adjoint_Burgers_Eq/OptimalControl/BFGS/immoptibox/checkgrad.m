function [maxJ, err, index] = checkgrad(fun, x, h, varargin)
%CHECKGRAD : Check the user's implementation of a nonlinear 
% (vector) function and its gradient (Jacobian).  
% The implementation must be a MATLAB function with declaration
%            function  [r, J] = fun(x,p1,p2,...)
% p1,p2,... are parameters of the function.  In connection with 
% nonlinear data fitting they may be arrays with coordinates of 
% the data points.   r  and  J  should be the vector function and 
% its Jacobian evaluated at  x .
%  
% Call  [maxJ, err, index] = checkgrad(fun, x, h, p1,p2,...)
%
% Input parameters 
% fun :  Handle to the function.
% x   :  The point where we wish to check the Jacobian.
% h   :  Steplength used in difference approximation to J .
%        May be a vector with  h(j)  used in the j'th coordinate
%        direction.
% p1,p2,...  are passed dirctly to the function FUN .
%
% Output parameters
% maxJ  :  The maximal absolute value of the elements in  J .
% err   :  Vector with three elements.  The maximal absolute deviation
%          between an element in  J  and the corresponding difference
%          approximation, computed by
%          err(1):  Forward difference with steplength  h
%          err(2):  Backward difference with steplength  h/2 
%          err(3):  "Extrapolated difference", see below. 
% index :  3*2 array, with  index(k,:)  giving the position in  J
%          where  err(k)  occurs. 
%
% If the function is "smooth" and  h  is sufficiently small (but not 
% so small that rounding errors dominate), then we can expect
%       error(1) "=" A*h,   error(2) "=" -A*(h/2)
% where "=" means "approximately equal to" and  A  is a constant, that
% depends on the function and  x , but not on  h .  The "extrapolated
% approximation" is  
%       (forward approximation + 2(backward approximation))/3
% and its error is of the order of magnitude  error(2)^2 .
%
% The algorithm is further discussed in "Checking Gradients", which
% can be downloaded from   http://www.imm.dtu.dk/~hbn/Software/#AUX

% Version 10.10.19.  hbn(a)imm.dtu.dk

%  Initialize
sx = size(x);   n = max(sx);
if  min(sx) ~= 1,  error('x  must be a vector'), end
sh = size(h);   nh = max(sh);
if  nh > 1 & (min(sh) ~= 1 | nh ~= n)
  error('h  must be a scalar or a vector of the same length as x')
end
if  nh == 1,  h = repmat(h,1,n); end
ieq = find(x(:)+h(:) == x(:) | x(:)-h(:)/2 == x(:));
if  ~isempty(ieq)
  error(['h was too small in direction(s)  ' int2str(ieq')]) 
end

% Check call
[f J] = feval(fun, x, varargin{:});                        % offset point
sf = size(f);   m = max(sf);
if  min(sf) ~= 1,  error('f  must be a vector'), end
sJ = size(J);  
if  m == 1  % Scalar case.  Allow both row and column
  ok = (min(sJ) == 1) & (max(sJ) == n);
  J = J(:).';
else        % Vector case.  Insist on  m*n
  ok = norm([m n] - sJ) == 0;
end
if  ~ok
  error('The sizes of  x , f  and  J  do not match'), end

% Check results
err = zeros(3,1);   index = zeros(3,2);
maxJ = max(max(abs(J)));   df = zeros(m,3);

% loop through the coordinate directions.
for j = 1 : n 
  %    Make differences for all function components
  xx = x;   xx(j) = x(j) + h(j);                   % forward point 
  dx = xx(j) - x(j);           % true difference in j'th direction        
  ff = feval(fun, xx, varargin{:});
  df(:,1) = (ff(:) - f(:))/dx;             % forward approximation
  xx(j) = x(j) - h(j)/2;   dx = x(j) - xx(j);     % backward point       
  fb = feval(fun, xx, varargin{:});
  df(:,2) = (f(:) - fb(:))/dx;            % backward approximation
  df(:,3) = (2*df(:,2) + df(:,1))/3;  % extrapolated approximation
  %  Store errors
  for  k = 1 : 3
    dif = df(:,k) - J(:,j); 
    [md, i] = max(abs(dif));    % largest difference and its index
    if  md > abs(err(k))                             % New maximum
      err(k) = dif(i);   index(k,:) = [i j];
    end
  end;
end  % j-loop
