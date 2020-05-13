%   IMMOPTIBOX  Optimization and Data Fitting Toolbox
%   Version 2.2, November 19, 2010
%   Copyright (c) 2010 by Hans Bruun Nielsen and IMM. 
% 
%   General unconstrained optimization
%     dampnewton - Damped Neton method.  Demands analytical expreesions
%                  for the gradient and Hessian
%     linesearch - Choice between soft and exact line search.
%     ucminf     - Unconstrained minimization of a scalar function of a
%                  vector variable.  Demands analytical expression for 
%                  the gradient.  Based on BFGS updating of the inverse                  
%                  Hessian and soft line search.
%
%   Unconstrained, nonlinear least squares problems
%     dogleg     - Powell's dog-leg method.  Demands analytical
%                  expression for the Jacobian.
%     marquardt  - Levenberg-Marquardt method.  Demands analytical
%                  expression for the Jacobian.
%     smarquardt - Levenberg-Marquardt method with successive updating
%                  of approximations to the Jacobian.   
% 
%   Data fitting with cubic splines
%     splinefit  - Weighted least squares fit of a cubic spline to
%                  given data points.  Possibility of assigning boundary
%                  conditions.  
%     splineval  - Evaluate a cubic spline s as computed by SPLINEFIT.
%     splinedif  - Evaluate s', s'' or s'''.
%
%   Robust estimation
%     huberobj    - Value and gradient of Huber estimator.
%                   Allows one-sided Huber function.
%     linhuber    - Minimizer of an extended linear Huber estimation
%                   problem.  Allows one-sided Huber function.
%     nonlinhuber - Minimizer of an extended nonlinear Huber estimation
%                   problem.  Allows one-sided Huber function.
%
%   Multiexponential fitting
%     mexpfit     - Weighted least squares fit of a multiexponential 
%                   model to given data points. Algorithm based on
%                   separability.
%
%   Nonlinear systems of equations
%     nonlinsys   - Solve nonlinear system of equations.  Dog Leg
%                   method with updating of approximate Jacobian.
%
%   Auxiliary programs
%     checkgrad   - Check implementation of gradient (or Jacobian)
%                   by means of finite differences.
% 
%   Test problems
%     uctpget    - Define test problem for unconstrained minimization.
%     uctpval    - Evaluate test problem.
%     
%   Data files
%     optic.dat  - Optic fibre data.
%     osl.dat    - Data from an optically stimulated luminiscence
%                  experiment.
%     peaks.dat  - Data with peaks and "shoulders".
%     wild.dat   - Data with "wild points".
%     efit1.dat  - Data for exponential fitting.
%     efit2.dat  - As efit1.dat, except that there are "wild points".