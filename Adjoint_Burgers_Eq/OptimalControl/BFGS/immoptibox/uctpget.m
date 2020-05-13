function [par, x0, tau0, delta0] = uctpget(p,m,n)
%UCTPGET  Define test problem for unconstrained minimization.
% Call:   [par, x0, tau0, delta0] = uctpget(p,m,n)
% Input parameters
% p   :  Problem number, integer in the range [1,22]
% m,n :  Number of components in the vectors f(x) and x,
%         respectively.  Not variable in all problems.
% Output parameters
% par    :  Struct defining the problem.
%           par.p  : Problem number 
%           par.pt = 0  signifies a least squares problem, (p <= 21)
%                    otherwise a general problem.
%           par.xm : Solution.  (NaN-s if m or n is free)
%           For some problems par has more fields.
% x0     :  Standard starting point.
% tau0   :  Standard value for the parameter opts(1) in the
%           Marquardt functions of the toolbox.
% delta0 :  Standard value for the parameter opts(1) in the
%           DOGLEG function of the toolbox.
%
% See H.B. Nielsen:  "UCTP - Test Problems for Unconstrained
% Optimization", Report IMM-REP-2000-18

% Version 04.01.31.  hbn(a)imm.dtu.dk

switch  p
  case 1    % Linear function.  Full rank
    par = struct('p',1, 'pt',0, 'A',eye(m,n) - (2/m)*ones(m,n), ...
      'xm',-ones(n,1));
    x0 = ones(1,n);   tau0 = 1e-8;   delta0 = 10;
  case 2    % Linear function.  Rank 1
    xm = (6/(2*m+1)/n/(n+1)) * ones(n,1);
    par = struct('p',2, 'pt',0, 'A',[1:m]'*[1:n], 'xm',xm);
    x0 = ones(1,n);   tau0 = 1e-8;   delta0 = 10;
  case 3    % Linear function.  Rank 1.  Zero cols. and rows
    xm = (6/(2*m-3)/(n-2)/(n+1)) * ones(n,1);
    par = struct('p',3, 'pt',0, 'A',[0 1:m-2 0]'*[0 2:n-1 0], 'xm',xm);
    x0 = ones(1,n);   tau0 = 1e-8;   delta0 = 10;
  case 4    % Rosenbrock
    par = struct('p',4, 'pt',0, 'xm',[1;1]); 
    x0 = [-1.2  1];  tau0 = 1;   delta0 = 1;
  case 5    % Helical Valley
    par = struct('p',5, 'pt',0, 'xm',[1;0;0]); 
    x0 = [-1  0  0];  tau0 = 1;   delta0 = 1;
  case 6    % Powell singular
    par = struct('p',6, 'pt',0, 'xm',zeros(4,1));
    x0 = [3  -1  0  1];   tau0 = 1e-8;   delta0 = 1;
  case 7    % Freudenstein and Roth
    xm = ([53;2] - sqrt(22)*[4;1])/3;
    par = struct('p',7, 'pt',0, 'xm',xm); 
    x0 = [.5  -2];  tau0 = 1;     delta0 = 1;
  case 8    % Bard
    Bard = [0.14  0.18  0.22  0.25  0.29  0.32  0.35  0.39  0.37 ... 
        0.58  0.73  0.96  1.34  2.10  4.39]';
    u = [1:15]';   v = 16 - u;   w = min([u'; v']).';
    xm = [8.241055996e-2;  1.133036099;  2.343695172];
    par = struct('p',8, 'pt',0, 'uvwy',[u v w Bard], 'xm',xm);
    x0 = ones(1,3);   tau0 = 1e-8;   delta0 = 1;
  case 9    % Kowalik and Osborne
    xm = [0.1928069346; 0.1912823287; 0.1230565069; 0.1360623307];
    y = [0.1957  0.1947  0.1735  0.1600  0.0844  0.0627 ...
        0.0456  0.0342  0.0323  0.0235  0.0246]';
    u = [4.0000  2.0000  1.0000  0.5000  0.2500  0.1670 ...
        0.1250  0.1000  0.0833  0.0714  0.0625]';
    par = struct('p',9, 'pt',0, 'yu',[y u], 'xm',xm);
    x0 = [.25  .39  .415  .39];  tau0 = 1;   delta0 = .1;
  case 10    % Meyer
    xm = [5.60963646990603e-003;  6.181346346e+003;  3.452236346e+002];
    par = struct('p',10, 'pt',0, 'ty',[45+5*(1:16)' Meyer], 'xm',xm);
    x0 = [.02  4e3  250];   tau0 = 1e-6;   delta0 = 100;
  case 11    % Watson
    t = [1:29]'/29;   B = ones(29,n);
    for  j = 2 : n,  B(:,j) = t.*B(:,j-1); end
    A = B(:,1:n-1) * diag(1:n-1);   xm = [];
    switch  n
      case  6,  xm = [-1.572508640e-2;  1.012434869; -2.329916260e-1;  
                       1.260430088; -1.513728923;  9.929964324e-1];
      case  9,  xm = [-1.530703649e-5;  9.997897039e-1;  1.476396369e-2
                        0.1463423283;  1.000821103; -2.617731141;
                        4.104403164;   -3.143612279;  1.052626408];      
      case 12,  xm = [-6.638060573e-9;  1.000001647; -5.639322310e-4
                       0.3478205409;   -0.1567315046;  1.052815176
                      -3.247271154;     7.288434898; -10.27184824
                       9.074113648;    -4.541375467;   1.012011889];
      otherwise, xm = repmat(NaN,n,1);
    end
    par = struct('p',11, 'pt',0, 'A',A, 'B',B, 'xm',xm);
    x0 = zeros(1,n);   tau0 = 1e-8;   delta0 = 1;
  case 12    % Box 3-d
    par = struct('p',12, 'pt',0, 't', .1*[1:m]', 'xm',[1; 10; 1]);
    x0 = [0  10  20];     tau0 = 1e-8;   delta0 = 1;
  case 13    % Jennrich and Sampson
    if  m == 10, xm = .25782521367*[1;1]; else, xm = [NaN; NaN]; end
    par = struct('p',13, 'pt',0, 't', [1:m]', 'xm',xm);
    x0 = [.3  .4];  tau0 = 1;   delta0 = .05;
  case 14    % Brown and Dennis
    if  m == 20
      xm = [-11.59443981; 13.20363001; -0.4034393057; 0.2367787344]; 
    else, xm = repmat(NaN,4,1); end
    par = struct('p',14, 'pt',0, 't', [1:m]'/5, 'xm',xm);
    x0 = [25  5  -5  -1];   tau0 = 1e-3;   delta0 = .5;
  case 15    % Chebyquad
    y = zeros(m,1);  
    if      (n == 8) & (m == 8)
      xm = [4.315283152e-2; 0.1930908899; 0.2663289000; 0.4999999000
             0.5000001000;  0.7336711000; 0.8069090000; 0.9568470101];
    elseif  (n == 8) & (m == 16)
      xm = [0.07974310852; 0.1783549801; 0.3101200148; 0.4412049499
             0.5587950200; 0.6898800200; 0.8216450649; 0.9202569150];
    elseif  (n == 9) & (m == 9)
      xm = [0.0442053010; 0.1994910050; 0.2356188900; 0.4160470950
        0.5000001100; 0.5839528850; 0.7643811500; 0.8005091500; 0.9557948250];
    elseif  (n == 9) & (m == 18)
      xm = [0.1018907415; 0.1894831408; 0.2944820487; 0.3970994868; 
        0.5000000000; 0.6029005132; 0.7055179513; 0.8105168592; 0.8981092585];
    else
      xm = repmat(NaN,n,1);
    end
    for  i = 2 : 2 : m,  y(i) = -1/(i^2 - 1); end
    par = struct('p',15, 'pt',0, 'y', y, 'xm',xm); 
    x0 = [1:n]/(n+1);  tau0 = 1;   delta0 = 1/(n+1);
  case 16    % Brown almost linear
    A = eye(n-1,n) + ones(n-1,n);
    par = struct('p',16, 'pt',0, 'A', A, 'xm',ones(n,1));
    x0 = .5*ones(1,n);  tau0 = 1;   delta0 = 1;
  case 17    % Osborne 1
    xm = [0.3754100521; 1.935846917;
      -1.464687141; 0.01286753465; 0.02212269965];
    par = struct('p',17, 'pt',0, 'ty',Osborne1, 'xm',xm);
    x0 = [.5  1.5  -1  .01  .02];     tau0 = 1e-8;   delta0 = .1;
  case 18    % Exponential fit.  n = 4
    xm = [-4.000026378; -4.999964833;  4.000242976; -4.000242571];
    par = struct('p',18, 'pt',0, 'ty',hbnefit, 'xm',xm);
    x0 = [-1 -2 1 -1];     tau0 = 1e-3;   delta0 = 1;
  case 19    % Exponential fit.  n = 2
    xm = [-4.000026708; -4.999964374];
    par = struct('p',19, 'pt',0, 'ty',hbnefit, 'xm',xm);
    x0 = [-1 -2];     tau0 = 1e-3;   delta0 = 1;
  case 20    % Scaled Meyer
    xm = [2.481778299; 6.181346346; 3.452236346];
    par = struct('p',20, 'pt',0, 'ty',[.45+.05*(1:16)' 1e-3*Meyer], 'xm',xm);
    x0 = [8.85 4 2.5];   tau0 = 1;   delta0 = 1;
  case 21    % Separated Meyer 
    xm = [6.181346324e+3;  3.452236339e+2];
    par = struct('p',21, 'pt',0, 'ty',[45+5*(1:16)' Meyer], 'xm',xm);
    x0 = [4e3  250];   tau0 = 1;   delta0 = 100;
   case 22    % Exp and squares 
     nn = [1:n].^2;   A = sum(1./nn);
     y = .75;   more = 1;
     while  more
       yo = y;   e = A*exp(-y);   y = y - (y-e)/(1+e);
       more = abs(y-yo) > 1e-12;
     end
     xm = exp(-y)./nn;
     par = struct('p',22, 'pt',1, 'xm',xm(:));
     x0 = zeros(1,n);   tau0 = 1e-3;   delta0 = 1;
  otherwise
    error(['Problem no ' int2str(p) '  is not available'])
  end
   
function  y = Meyer
  y = [34780  28610  23650  19630  16370  13720  11540 ...
        9744   8261   7030   6005   5147   4427   3820 ...
        3307   2872]';
       
function  ty = Osborne1
  y = [0.844  0.908  0.932  0.936  0.925  0.908  0.881 ...
       0.850  0.818  0.784  0.751  0.718  0.685  0.658 ...
       0.628  0.603  0.580  0.558  0.538  0.522  0.506 ...
       0.490  0.478  0.467  0.457  0.448  0.438  0.431 ...
       0.424  0.420  0.414  0.411  0.406]';
  t = 10 * (0 : 32);  ty = [t(:) y];
       
function  ty = hbnefit
  y =[0.090542  0.124569  0.179367  0.195654  0.269707 ...
      0.286027  0.289892  0.317475  0.308191  0.336995 ...
      0.348371  0.321337  0.299423  0.338972  0.304763 ...
      0.288903  0.300820  0.303974  0.283987  0.262078 ...
      0.281593  0.267531  0.218926  0.225572  0.200594 ...
      0.197375  0.182440  0.183892  0.152285  0.174028 ...
      0.150874  0.126220  0.126266  0.106384  0.118923 ...
      0.091868  0.128926  0.119273  0.115997  0.105831 ...
      0.075261  0.068387  0.090823  0.085205  0.067203]';
  t = .02 * [1:45];   ty = [t(:) y];
