function [u, u1] = minNewtonCG(u,y0,z,Opts)

global tol_line min_alpha eta N Nt time space

max_opt       = Opts.max_opt;        % optimization iterations     
tol_opt_f     = Opts.tol_opt_f;      % change in objective function 
tol_opt_gradf = Opts.tol_opt_gradf;  % tolerance for gradient
tol_line = Opts.tol_line;            % tolerance for line search
min_alpha = Opts.min_alpha;          % minimum step size in direction s_k
eta = Opts.eta;                      % tolerance for truncated CG


U = reshape(u,N-1,Nt);
Yfull  = Burgers(U,y0);

for k = 1:max_opt
    tic;
    
    cost_old = costFun(u,Yfull,z);

    %% Compute grad_f via adjoint I 
    [Lambda, gradf] = Adjoint_Grad(Yfull, U, z);

    % Check for stopping of optimization
    if (norm(gradf) < tol_opt_gradf) % first order optimality condition
        disp('Zero-gradient criterium satisfied.');
        return
    end    

    %% Solve Newton eqn. for s_k via CG alg, Hess_f*v is calculated via adjoint II
    s0 = zeros((N-1)*Nt,1);
    [s,iter] = TruncatedCG( @(u)Adjoint_HessVec(u,Yfull,Lambda), s0, -gradf, (N-1)*Nt );

    %% Perform line search in s_k direction
    [alpha, Yfull, u] = ArmijoLineSearch( @costFun, Yfull, y0, z, u, s, s'*gradf );    
    if k == 1
        u1 = u;
    end
    U = reshape(u,N-1,Nt); % update

    cost_new=costFun(u,Yfull,z); 

    if (abs(cost_new-cost_old) < tol_opt_f)
        disp('No further improvement in objective function');
        return
    end
    
    count = toc;
    PlotAndOutput(time, space, Yfull, z, Lambda, U, k, gradf, s, alpha, count, iter);
    
end




end

