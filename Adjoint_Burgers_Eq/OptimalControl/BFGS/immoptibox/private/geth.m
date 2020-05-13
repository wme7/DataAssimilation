function  [h, mu] = geth(A,g,mu)
% Solve  (A + mu*I)h = -g  with possible adjustment of  mu

% Version 10.02.17.  hbn@imm.dtu.dk

mA = max(abs(A(:)));
if  mA == 0  % zero matrix
  h = g / (-mu);  return
end

% Factorize with check of pos. def.
n = size(A,1);  chp = 1;
while  chp
  [R chp] = chol(A + mu*eye(n));
  if  chp == 0  % check for near singularity
    chp = rcond(R) < 1e-15;
  end
  if  chp,  mu = max(10*mu, eps*mA); end
end

% Solve  (R'*R)h = -g
h = R \ (R' \ (-g));   