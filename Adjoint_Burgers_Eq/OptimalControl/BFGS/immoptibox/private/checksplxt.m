function  [err, tyw, ma] = checksplxt(A,B,cs,caller)
%CHECKSPLXT Check spline knots and data abscissae.

% Version 04.05.05.  hbn@imm.dtu.dk

err = 0;   tyw = [];  ma = 0;
if  cs == 1  % Check knots and data
  x = A;  % knots
  sx = size(x);   
  if  min(sx) ~= 1 | ~isreal(x) | any(isnan(x(:)) | isinf(norm(x(:)))) | max(sx) < 2
    err = 1; 
  else
    d = diff(x);
    if  any(d < 0),  err = 2;
    elseif  d(1) == 0 | d(end) == 0,  err = 6; end
  end
  if  ~err  % Check and order data
    [m q] = size(B);
    if      m == 0,  err = 5;
    elseif  q < 2,   err = 4;
    else
      t = B(:,1);
      if  ~isreal(t) | any(isnan(t) | isinf(t) | t < x(1) | t > x(end))
        err = 3; 
      else  % order data
        if  q > 2
          j = find(~isinf(B(:,3)) & ~isnan(B(:,3)) & ...
            imag(B(:,3)) == 0 & B(:,3) > 0);   
          ma = length(j);
          if  ma < m,  tyw = B(j,1:3); 
          else,        tyw = B(:,1:3); end
        else
          ma = m;   tyw = [B ones(m,1)]; 
        end   % compress data
        
        % sort data
        [t js] = sort(tyw(:,1));  
        tyw = tyw(js,:);
      end      
    end % data ok
  end % knots ok
  
else  % A is a spline and B holds arguments
  if  A.fail ~= 0
    error([caller ':  The spline is not defined']),  end
  if  min(size(B)) ~= 1,  err = 1; end
  if  ~err  % check vector elements
    if  ~isreal(B) | any(isnan(B) | isinf(B)),  err = 1; end
  end
  if  err
    error([caller ':  t is not a real valued vector']), end
  
  if  cs > 2  % Check that t is in knot range
    if  any(B < A.x(1) | B > A.x(end))
      error([caller ':  t is not in the knot range']), end
  end
end