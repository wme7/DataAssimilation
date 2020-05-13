function [xk,info] =  spg_cons(func,Opts)

global flag_fg

% Matlab version of SPG:  Birgin, Martinez and Raydan (SIAM J. Opt. 2000)
%
% Possible ot use of parameter eta in LINE SEARCH as in
% Martinez, La Cruz and Raydan, (Math. of Comput. 2006)
%
% The user must supply func and grad for the function and the gradient.
%
% Author:      M. Raydan et al.
%
% Modified by: M. Rojas
% Date:       15 Mar 2011, changed interface to allow different functions
% Date:       16 Mar 2011, added 'gamma' and 'alpha' as fields of 'Opts'
% Date:       27 Nov 2011, changed fmax to an m-vector, fixed bug related to m
% Date:        9 Apr 2013, unconstrained version
%


% spg setup - begin

% func = F.func;
% if isfield(F,'param')
%    funcparam = F.param;
% else
%    funcparam = [];
% end
% grad      = F.grad; % for Slawek & Manuel: func returns f&g at x


n = Opts.n;
if isfield(Opts,'Proj')
   proj = Opts.Proj;
else
   proj = @(x)x; %identity as default
end
if isfield(Opts,'x0')
   xk = Opts.x0;
else
   xk = ones(n,1);
end
xk = proj(xk);

if isfield(Opts,'tol')
   tol = Opts.tol;
else
   tol = 1e-6;
end
if isfield(Opts,'maxiter')
   maxiter = Opts.maxiter;
else
   maxiter = 200;
end

if isfield(Opts,'rhomin')
   rhomin = Opts.rhomin;
else
   rhomin = 1e-15;
end
if isfield(Opts,'rhomax')
   rhomax = Opts.rhomax;
else
   rhomax = 1e+15;
end
if isfield(Opts,'res')
   res = Opts.res;
else
   res = 1000;
end
if isfield(Opts,'maxfc')
   maxfc = Opts.maxfc;
else
   maxfc = 3000;
end
if ~isfield(Opts,'gamma')
   Opts.gamma = 1e-4;
end
if ~isfield(Opts,'alpha')
   Opts.alpha = 1.0;
end
if isfield(Opts,'m')
   m = Opts.m;
else
   m = 2;
end
if isfield(Opts,'xtrue')
   xtrue = Opts.xtrue;
else
   xtrue = ones(n,1);
end

if isfield(Opts,'err_f')
    err_f = Opts.err_f;
else
    err_f  = 1e-8;
end
if isfield(Opts,'err_x')
    err_x = Opts.err_f;
else
    err_x  = 1e-8;
end

% spg setup - end

% Two key Parameters
% m     = Nonmonotone parameter in LS  
% eta   = Initial value of parameter for additional nonmonotone behavior 

    iter  = 0;  
    fcnt  = 0;
    fxk   = ones(m,1);
    gcnt  = 0;
    backt = 0;

    flag_fg = 1;
    [f,g]    = func(xk);
    fcnt = fcnt + 1;
    gcnt = gcnt + 1;
    
    fxk      = f*fxk;
    
    if isfield(Opts,'eta')
       if ischar(Opts.eta) & strcmp(Opts.eta,'f') 
          Opts.eta  = abs(f);
       end
    else
       Opts.eta  = 0;
    end

    % COMPUTE CONTINUOUS PROJECTED GRADIENT (AND ITS NORM)
    pg      = xk-g;
    pg      = proj(pg);
    pg      = pg-xk;
    
    pgnorm  = norm(pg,'inf');

    % DEFINE INITIAL SPECTRAL STEPLENGTH
    if (pgnorm ~= 0)
       rho = min(rhomax,max(rhomin,1.0d0/pgnorm));
    end
         
%     err_f  = 1e-8;
%     err_x  = 1e-8;
    f_norm = err_f ;
    x_norm = err_x ;
    fprintf('  iter | fcnt | gcnt| f_norm | x_norm | restemp |        f       |\n');
       
    while (pgnorm >= tol)   & ...
          (iter <= maxiter) & ...
          (fcnt <= maxfc)   & ...
          (f_norm >= err_f) & ...
          (x_norm >= err_x)

       iter = iter + 1;
        relspgk(iter) = norm(xtrue-xk)/norm(xtrue);
        itnorm(iter)  = norm(xk);
        fback = f;
       d    = xk - rho*g;
       d    = proj(d);
       d    = d - xk;
       gtd  = g'*d;

       % CALL  LINE SEARCH   
%%       [fcnt,fnew,xnew,back,fxk] = ...
%%        linesearch(iter,xk,gtd,fcnt,maxfc,m,d,f,fxk,eta,func,funcparam);
      
       [fnew,xnew,back] = linesearch;
       flag_fg = 3;
       [~,gnew] = func(xnew);
       gcnt = gcnt + 1;
       
       f         = fnew;
       indx      = mod(iter,m)+1;
       fxk(indx) = f;
       backt     = backt + back;

    % COMPUTE S = XNEW-X, Y = GNEW-G, <S,S>, <S,Y>
       S       = xnew - xk;
       Y       = gnew - g;
       sts     = S'*S;
       sty     = S'*Y;
       xk      = xnew;

       %x_k_all(1:300,iter) = xnew;
       g       = gnew;
       pg      = xk - g;
       pg      = proj(pg); 
       pgnorm  = norm(pg-xk,'inf');
       if (sty <= 0)
           rho = rhomax;
       else
           rho = min(rhomax,max(rhomin,sts/sty));
       end
    
       %restemp = norm(pg-xk,'inf');
        restemp = norm(pg-xk);

      f_norm = norm(f-fback);
      x_norm =  norm(S);
       
%************************************************
      res = [res restemp];

      fprintf('  %5d | %5d | %5d | %4.2e | %4.2e  | %4.2e | %5f \n',iter,fcnt,gcnt,f_norm,x_norm,restemp, f);
      
%************************************************    
    
    end
   

   
   relspgk(iter) = norm(xtrue-xk)/norm(xtrue);
   itnorm(iter)  = norm(xk);

   info.iter    = iter;
   info.f       = f;
   info.fcnt    = fcnt;
   info.gcnt    = gcnt;
   info.backt   = backt;
   info.pgnorm  = pgnorm;
   info.relspgk = relspgk;
   info.itnorm  = itnorm;
   info.stop    = [pgnorm < tol; ...
                   iter > maxiter; ...
                   fcnt > maxfc; ...
                   f_norm < err_f ; ...
                   x_norm < err_x];
   info.res     = res;


%************************************************
%     Nested Functions 
%************************************************    
function [fnew,xnew,back] = linesearch

   gamma = Opts.gamma;
   alpha = Opts.alpha;
   eta   = Opts.eta;

   back = 0;
   etak = eta/(iter^1.1);

   fmax = max(fxk);

   xnew = xk + d;
   
   flag_fg = 2;
   [fnew,~] = func(xnew);
   fcnt        = fcnt + 1; 
   
   while (fnew > fmax + (gamma*alpha*gtd) + etak) & (fcnt < maxfc),
      if (alpha <= 0.1)
         alpha = alpha/2;
      else
         alphatemp = (-gtd*alpha^2)/(2*(fnew - f - (alpha*gtd)));
         if (alphatemp < 0.1) | (alphatemp > 0.9*alpha)
            alphatemp = alpha/2;
         end
         alpha = alphatemp;
      end
      xnew = xk + (alpha*d);
      
      flag_fg = 2;
      [fnew,~] = func(xnew);
      fcnt = fcnt + 1;

      back = back + 1;
      
   end
   %gnew = grad_test(xnew);

end % linesearch
end % spg


