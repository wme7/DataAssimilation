function  [F, dF, d2F] = uctpval(x,par)
%UCTPFGH  Evaluate test problem for nconstrained minimization,
% as defined by UCTPGET.
% Call
%    F = uctpval(x,par)
%    [F, dF] = uctpval(x,par)
%    [F, dF, d2F] = uctpval(x,par)
% Input parameters
% x   :  Argument.
% par :  Struct defining the problem.
% Output parameters
% F   :  The function.  If  par.p <= 21 and par.pt = 0,
%        then  F  is a vector, otherwise a scalar.
% dF  :  If  par.p <= 21 and par.pt = 0, then  dF  is the
%        Jacobian J(x), otherwize the gradient  F'(x) .
% d2F :  Presumes par.pt ~= 0.  Hessian matrix

% Version 04.04.15.  hbn(a)imm.dtu.dk

% Check for Hessian
if  nargout > 2,  par.pt == 1; end

% Ensure that  x  is a column vector
x = x(:);

if  par.p <= 21    % Least squares problem
  difH = 0;  % Not difference approximation to Hessian
  switch  par.p
    case {1,2,3}    % Linear function
      J = par.A;  f = J*x - 1;  
      if  nargout > 2,   H = J'*J; end
    case 4    % Rosenbrock
      f = [10*(x(2) - x(1)^2); 1-x(1)]; 
      if  nargout > 1
        J = [-20*x(1)  10; -1  0];  
        if  nargout > 2,  H = J'*J + f(1)*[-20 0; 0 0]; end
      end
    case 5    % Helical Valley
      t = atan(x(2)/x(1))/(2*pi);
      if  x(1) < 0,  t = t + .5; end
      nx = norm(x(1:2));  nx2 = nx*nx;
      f = [10*(x(3) - 10*t); 10*(nx - 1); x(3)];
      if  nargout > 1
        K1 = 50/pi/nx2;   K2 = 10/nx;
        J = [K1*x(2)  -K1*x(1)  10; K2*x(1:2)' 0; 0 0 1];
        if  nargout > 2
          H = J'*J;
          q1 = x(1)^2;  q2 = x(2)^2;  p = x(1)*x(2);  ii = 1:2;
          H(ii,ii) = H(ii,ii) + f(1)*K1/nx2*[-2*p q1-q2; q1-q2 2*p] ...
            + f(2)*K2/nx2*[q2 -p; -p q1];
        end
      end
    case 6    % Powell singular
      s5 = sqrt(5);   s10 = sqrt(10);
      d3 = x(2) - 2*x(3);   d4 = x(1) - x(4);
      f = [x(1)+10*x(2); s5*(x(3) - x(4)); d3^2; s10*d4^2];
      if  nargout > 1
        J = [1 10 0 0; 0 0 s5 -s5; 0 2*d3*[1 -2] 0; 2*s10*d4*[1 0 0 -1]];
        if  nargout > 2
          H = J'*J;
          H(2:3,2:3) = H(2:3,2:3) + f(3)*[2 -4;-4 8];
          H([1 4],[1 4]) = H([1 4],[1 4]) + f(4)*2*s10*[1 -1;-1 1];
        end
      end
    case 7    % Freudenstein and Roth
      x1 = x(1);   x2 = x(2);
      f = [(x1 - x2*(2 - x2*(5 - x2)) - 13)
        (x1 - x2*(14 - x2*(1 + x2)) - 29)];
      if  nargout > 1
        J = [1  (-2 + x2*(10 - 3*x2)); 1  (-14 + x2*(2 + 3*x2))];
        if  nargout > 2
          H = J'*J;
          H(2,2) = H(2,2) + f(1)*10-6*x2 + f(2)*(2+6*x2);
        end
      end
    case 8    % Bard
      D = par.uvwy(:,2:3)*x(2:3);
      f = par.uvwy(:,4) - (x(1) + par.uvwy(:,1)./D);
      if  nargout > 1
        F = -par.uvwy(:,1)./D.^2;
        J = -[ones(size(D))  (F*[1 1]).*par.uvwy(:,2:3)];
        if  nargout > 2
          H = J'*J;   ii = 2:3;   A = zeros(2,2);
          for  i = 1 : length(f)
            u = par.uvwy(i,1);  v = par.uvwy(i,2);  w = par.uvwy(i,3);
            J(i,:) = [-1 [v w]*u/D(i)^2];
            A = A + f(i)*[-2*v^2  v*w; v*w -2*w^2]*u/D(i)^3;
          end
          H(ii,ii) = H(ii,ii) + A;
        end
      end
    case 9    % Kowalik and Osborne
      u = par.yu(:,2);   x1 = x(1);
      N = u.*(u + x(2));   D = u.*(u + x(3)) + x(4);
      f = par.yu(:,1) - x1*N./D;
      if  nargout > 1
        F = -x1*N./D.^2;
        J = -[N./D  x1*u./D  F.*u  F];
        if  nargout > 2
          H = J'*J;
          for  i = 1 : 11
            Hi = zeros(4,4);
            T = x(1)*N(i);  ui = u(i);  Di = D(i);  Ni = N(i);  u1 = [ui 1];
            a = [-1 u1*(Ni/Di)]*(ui/Di);
            Hi(1,2:4) = a;  Hi(2:4,1) = a';
            Hi(2,3:4) = u1*(x1/Di^2);  Hi(3:4,2) = Hi(ii(2),3:4)';
            a = -2*T/Di^3;   
            Hi(3,3:4) = u1*(ui*a);  Hi(4,3) = Hi(ii(3),4);
            Hi(4,4) = a;   H = H + f(i)*Hi;
          end
        end
      end
    case 10    % Meyer
      D = par.ty(:,1) + x(3);   q = x(2)./D;  e = exp(q);
      f = x(1)*e - par.ty(:,2);
      if  nargout > 1
        J = [e  x(1)*e./D  -x(1)*x(2)*e./D.^2];
        if  nargout > 2
          H = J'*J;          
          for  i = 1 : m
            p = x(1)/D(i);  a = q(i);  b = p*a;  c = -p*(1+a);
            H = H + f(i)*(e(i)/D(i))*[0 1 -a;1 p c;-a c b*(2+a)];
          end
        end
      end
    case 11    % Watson
      m = 31;  n = length(x);
      A = [zeros(29,1) par.A];   B = par.B;
      g = B*x;   x1 = x(1);   x2 = x(2);  
      f = [(par.A*x(2:end) - g.^2 - 1); x1; x2-x1^2-1];
      if  nargout > 1
        J = zeros(31,n);   g = -2*g;
        J(:,1) = [g; 1; -2*x1];   J(31,2) = 1;
        for  j = 2 : n
          J(1:29,j) = par.A(:,j-1) + g.*par.B(:,j);
        end
        if  nargout > 2
          H = J'*J;
          for  i = 1 : 29
            bi = B(i,:);
            H = H - (2*f(i)*bi')*bi;  
          end
          H(1,1) = H(1,1) - 2*f(end);
        end
      end
    case 12    % Box 3-d
      t = par.t;   E = exp(-t*[x(1:2)' 1  10]);
      c = [1; -1; -x(3)*[1; -1]];
      f = E*c;
      if  nargout > 1
        J = [-t.*E(:,1) t.*E(:,2) E(:,3:4)*[-1;1]];
        if  nargout > 2
          H = J'*J;  ii = 1:2;
          for  i = 1 : length(t)
            H(ii,ii) = H(ii,ii) + f(i)*t(i)*diag(J(i,ii));  
          end
        end
      end
    case 13    % Jennrich and Sampson
      t = par.t;   E = exp(t*x');
      f = 2*(1 + t) - E*[1; 1];
      if  nargout > 1
        J = -(t*ones(1,2)).*E;
        if  nargout > 2
          H = J'*J;  
          for  i = 1 : length(t)
            H = H + f(i)*t(i)*diag(J(i,:));  
          end
        end
      end
    case 14    % Brown and Dennis
      t = par.t;   st = sin(t);
      d1 = x(1) + x(2)*t - exp(t);   d2 = x(3) + x(4)*st - cos(t);
      f = d1.^2 + d2.^2;
      if  nargout > 1
        J = 2*[d1  d1.*t  d2  d2.*st];
        if  nargout > 2
          H = J'*J;  
          for  i = 1 : length(t)
            a = [1 t(i) 1 st(i)];
            H = H + a'*(2*f(i)*a);
          end
        end
      end
    case 15    % Chebyquad
      z = real(acos(2*x' - 1));   f = -par.y;
      m = length(f);   n = length(x);
      for  r = 1 : m,  f(r) = sum(cos(r*z))/n + f(r); end
      if  nargout > 1
        J = zeros(m,n);   d = sqrt(1 - (2*x' - 1).^2);
        nz = find(d);   iz = find(d == 0);
        for  r = 1 : m
          J(r,nz) = 2/n*r*sin(r*z(nz)) ./ d(nz);
          J(r,iz) = 2/n*r^2;
        end 
        J = real(J);
        if  nargout > 2,  difH = 1; end
      end
    case 16    % Brown almost linear
      n = length(x);   p = prod(x);
      f = [par.A*x-(n+1); p-1];
      if  nargout > 1
        pp = zeros(1,n);
        for  j = 1 : n,  pp(j) = prod(x([1:j-1 j+1:n])); end 
        J = [par.A; pp];
        if  nargout > 2  % Hessian
          if      n == 1,  H = 0;
          elseif  n == 2,  H = [0 1;1 0];
          else
            for  i = 1 : n-1
              for  j = i+1 : n
                H(i,j) = prod(x([1:i-1 i+1:j-1 j+1:n]));
                H(j,i) = H(i,j);
              end
            end
          end
        end
      end
    case 17    % Osborne 1
      t = par.ty(:,1);   E = exp(-t*x(4:5)');
      f = par.ty(:,2) - x(1) - E*x(2:3);
      m = length(t);  n = length(x);
      if  nargout > 1
        J = -ones(length(t),5);
        J(:,2:5) = -[E  -x(2)*t.*E(:,1)  -x(3)*t.*E(:,2)];
        if  nargout > 2
          A4 = zeros(5,5);
          A1 = A4;  A1(2,4) = 1;  A1(4,2) = 1;
          A2 = A4;  A2(3,5) = 1;  A2(5,3) = 1;
          A3 = A4;  A3(4,4) = -x(2);  A4(5,5) = -x(3);
          H = J'*J;
          for  i = 1 : length(t)
            H = H + f(i)*t(i)*(E(i,1)*(A1 + t(i)*A3) ...
              + E(i,2)*(A2 + t(i)*A4));
          end
        end
      end
    case 18    % Exponential fit.  n = 4
      t = par.ty(:,1);   E = exp(t*x(1:2)');
      f = par.ty(:,2) - E*x(3:4);
      if  nargout > 1
        J = -[E  E];
        J(:,1:2) = -[x(3)*t.*E(:,1)   x(4)*t.*E(:,2)];
        if  nargout > 2
          H = J'*J;  
          j1 = [1 3];  j2 = [2 4];  x3 = -x(3);  x4 = -x(4);
          for  i = 1 : length(t)
            H(1,j1) = H(1,j1) + f(i)*t(i)*E(i,1)*[x3*t(i) 1];
            H(2,j2) = H(2,j2) + f(i)*t(i)*E(i,2)*[x4*t(i) 1];
          end
          H(3,1) = H(1,3);  H(4,2) = H(2,4);
        end
      end
    case 19    % Exponential fit.  n = 2
      t = par.ty(:,1);   E = exp(t*x');   c = E\par.ty(:,2);
      f = par.ty(:,2) - E*c;
      if  nargout > 1
        A = E'*E;   H = [t t] .* E;
        G = A\ (diag(H'*f) - H'*E*diag(c));
        J = -(E*G + H*diag(c));
        if  nargout > 2,  difH = 1; end
      end
    case 20    % Modified Meyer
      D = par.ty(:,1) + x(3);   e = exp(10*x(2)./D - 13);
      f = x(1)*e - par.ty(:,2);
      if  nargout > 1
        q = (10*x(1)./D) .* e;
        J = [e  q  -(x(2)./D).*q];
        if  nargout > 2
          H = J'*J;  A = zeros(3,3);   
          x1 = x(1);  x2 = x(2);
          for  i = 1 : length(f)
            Di = D(i);   A(1,:) = [0 Di -1];
            A(2,:) = [Di 10*x1^2 -x1*(1 +10*x2/Di)];
            A(3,:) = [-1 -x1*(1 +10*x2/Di) 2*x1*(1 +5*x2/Di)];
            H = H + f(i)*e(i)/Di^2*A;
          end
        end            
      end
    case 21    % Separated Meyer
      x2 = par.ty(:,1) + x(2);   a = x(1) ./ x2;
      F = exp(a);   E = F'*F;   c = (F' * par.ty(:,2))/E;
      f = F*c - par.ty(:,2);
      if  nargout > 1
        d1 = x(1) ./ x2;   d2 = c ./ x2;
        by = (f + F*c) ./ x2;   by = [by -by.*d1];   
        dc = -(F'*by)/E;
        J = [F F] .* [dc(1)+d2  dc(2)-d1.*d2];
        if  nargout > 2,  difH = 1; end
      end
      
  end
    
  % Reformulate to scalar problem
  if  par.pt == 0  % Vector function
    F = f;  
    if  nargout > 1,  dF = J; end  % Jacobian
  else  % Scalar function
    F = .5 * norm(f)^2;
    if  nargout > 1  % Gradient 
      dF = J'*f;
      if  nargout > 2  % Hessian
        if  difH  % difference approximation
          n = length(x);   H = zeros(n,n);   u = x;
          for  j = 1 : n
            d = 1e-5*(abs(x(j)) + 1e-5);  % step length
            fj = x(j) + d;   bj = x(j) - d;
            u(j) = fj;  [dum gf] = uctpval(u,par);
            u(j) = bj;  [dum gb] = uctpval(u,par);
            H(:,j) = (gf - gb)/(fj - bj);
          end
        end
        H = (H + H')/2;  % symmetrize
        d2F = H;
      end
    end
  end
  
else    % 'Born' scalar function
  switch  par.p
    case 22    % Exp and squares
      n = length(x);
      e = exp(-sum(x));  ii = (1:n).^2;
      F = e + .5*sum(ii(:) .* x.^2);
      if  nargout > 1 
        dF = -e + ii(:) .* x; 
        if  nargout > 2
          d2F = diag(ii) + e;
        end
      end
  end
end
