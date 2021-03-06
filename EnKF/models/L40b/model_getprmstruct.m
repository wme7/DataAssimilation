% function [prm, fmt] = model_getprmstruct
%
% Reads custom parameter structure for L40p model.
%
% @return prm - system parameters
% @return fmt - parameter formats
%
% Description of individual fields in the parameter file:
%
% field = approx
%   - flag; whether to use native Matlab ODE integrator with step control (if
%     0) or make one integrator step per one model step (if 1)
% values =
%   0
%   1*
%
% field = integrator
%   - integrator tag
% values =
%   rk4 - for Runge-Kutta 4-order integrator
%   dp5 - for Dorman-Prince 5-order integrator
%
% field = dt
%   - time step for one model step
% values =
%   <value> (0.05*)
%
% field = F
%   - value of forcing parameter in L40 model
% values =
%   <value> (8.0*)
%
% field = T
%   - period (in model steps) of change of bias
% values =
%   <value> (50)
%
% field = use_compiled
%   - flag; whether ot use compiled C code or Matlab code for stepping the
%     model
% values =
%   0
%   1*
%
% field = fname_samples
%   - name of the file with samples for initial ensemble
% values =
%   <file name>

% File:           model_getprmstruct.m
%
% Created:        25/08/2010
%
% Last modified:  25/08/2010
%
% Author:         Pavel Sakov
%                 NERSC
%
% Last modified:  25/08/2010
%
% Purpose:        Creates two structures: prm -- to hold the system parameters;
%                 and fmt -- to hold the format strings for each parameter,
%                 used while reading from a file
%
% Description:    
%
% Revisions:

%% Copyright (C) 2010 Pavel Sakov
%% 
%% This file is part of EnKF-Matlab. EnKF-Matlab is a free software. See 
%% LICENSE for details.

function [prm, fmt] = model_getprmstruct
    
    prm = struct(                            ...
        'integrator',    {'dp5'},            ...
        'dt',            {0.05},             ...
        'T',             {50.0},             ...
        'F',             {8},                ...
        'use_compiled',  {1},                ...
        'fname_samples', {'L40_samples.mat'} ...
        );
        
    fmt = struct(                            ...
        'integrator',    {'%s'},             ...
        'dt',            {'%g'},             ...
        'T',             {'%g'},             ...
        'F',             {'%g'},             ...
        'use_compiled',  {'%d'},             ...
        'fname_samples', {'%s'}              ...
        );

    return
