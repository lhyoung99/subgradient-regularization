function options = setdefaults_gradsamp(dimension, options)
%  to be called only by gradsamp
%  check that fields of pars and options are set correctly and
%  set basic default values for options for gradsamp
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  HANSO 3.0 Copyright (C) 2021  Michael Overton
%%  This program is free software: you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation, either version 3 of the License, or
%%  (at your option) any later version.
%%
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%
%%  You should have received a copy of the GNU General Public License
%%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% fields of pars
%
% if ~isfield(pars, 'nvar')
%    error('setdefaults_gradsamp: input "pars" must have field "nvar"')
% elseif ~isposint(pars.nvar)
%    error('setdefaults_gradsamp: input "pars.nvar" (number of variables) must be a positive integer')
% end
% if ~isfield(pars,'fgname')
%    error('setdefaults_gradsamp: input "pars" must have field "fgname" (name of m-file computing function and gradient)')
% end
% 
% fields of options
%
if isfield(options, 'prtlevel')
    if options.prtlevel ~= 0 && options.prtlevel ~= 1 && options.prtlevel ~= 2 &&...
            options.prtlevel ~= 3
        error('setdefaults_gradsamp: input "options.prtlevel" must be 0,1,2, or 3')
    end
else
    options.prtlevel = 1;
end

if isfield(options, 'print_frequency')
    if ~isposint(options.print_frequency)
        error('setdefaults_gradsamp: input "options.print_frequency" must be positive integer')
    end
else
    options.print_frequency = 10000;
end

% if isfield(options,'x0')
%     if any(size(options.x0) ~= [pars.nvar 1])
%         error('setdefaults_gradsamp: input "options" must have field "x0" with size "pars.nvar" by 1')
%     end
% else
%     fprintf('setdefaults_gradsamp: setting options.x0 randomly\n')
%     options.x0 = randn(pars.nvar,1);
% end
if isfield(options,'eps0') % initial value of sampling radius
    if ~isposreal(options.eps0)||isinf(options.eps0)
        error('setdefaults_gradsamp: input "options.eps0" must be positive')
    end
else
    options.eps0 = 0.1;
end
if isfield(options,'nu0') % initial stationarity target
    if ~isposreal(options.nu0)
        error('setdefaults_gradsamp: input "options.nu0" must be positive')
    end
else
    options.nu0 = 0.1;
end
if isfield(options,'ngrad')
    if ~isnonnegint(options.ngrad) % could be 0 to simulate steepest descent
        error('setdefaults_gradsamp: input "options.ngrad" must be a nonnegative integer')
    end
    if options.ngrad == 0 && options.prtlevel > 0
        fprintf('setdefaults_gradsamp: warning: number of sampled gradients is zero: steepest descent\n')
    elseif options.ngrad <= dimension + 1 && options.prtlevel > 0
        fprintf('setdefaults_gradsamp: warning: number of sampled gradients < number of variables + 1\n')
    end
else
    options.ngrad = 2*dimension;
end
if isfield(options,'beta') % line search Armijo-type condition
    if ~isposreal(options.beta) || options.beta >= 1
        error('setdefaults_gradsamp: input options.beta must be in (0,1)')
    end
else
    options.beta = 1e-4;
end
if isfield(options,'gamma') % line search contraction parameter
    if ~isposreal(options.gamma) || options.gamma >= 1
        error('setdefaults_gradsamp: input options.gamma must be in (0,1)')
    end
else
    options.gamma = 0.5;
end
if isfield(options,'eps_opt') % sampling radius termination tolerance
    if ~isposreal(options.eps_opt)
        % Alg GS actually allows eps_opt = 0 but that's only for
        % theoretical convergence analysis, it makes no sense in
        % the implementation
        error('setdefaults_gradsamp: input options.eps_opt must be positive')
    end
else
    options.eps_opt = 1e-6;
end
if options.eps0 < options.eps_opt
    error('setdefaults_gradsamp: input options.eps0 must be >= options.eps_opt')
end
if isfield(options,'nu_opt') % stationarity termination tolerance
    if ~isnonnegreal(options.nu_opt)
        % Alg GS actually allows nu_opt = 0 but that's only for
        % theoretical convergence analysis, it makes no sense in
        % the implementation
        % [Modified in 2025:] allow nu_opt = 0 since
        % a new stopping condition on obj_tol has been included
        error('setdefaults_gradsamp: input options.nu_opt must be nonnegative')
    end
else
    options.nu_opt = 1e-6;
end
if options.nu0 < options.nu_opt
    error('setdefaults_gradsamp: input options.nu0 must be >= options.nu_opt')
end
if isfield(options,'theta_eps') % reduction factor for sampling radius
    if ~isposreal(options.theta_eps) || options.theta_eps > 1
        error('setdefaults_gradsamp: input options.theta_eps must be in (0,1]')
    end
else
    options.theta_eps = 0.1;
end
if isfield(options,'theta_nu') % reduction factor for stationarity target
    if ~isposreal(options.theta_nu) || options.theta_nu > 1
        error('setdefaults_gradsamp: input options.theta_nu must be in (0,1]')
    end
else
    options.theta_nu = 0.1;
end
if isfield(options, 'delta') % termination parameter for line search
                             % not in algorithm statement, but needed
                             % because of rounding errors
    if ~isposreal(options.delta)
        error('setdefaults_gradsamp: input "options.delta" must be positive')
    end
else
    options.delta = 1e-16;
end
if isfield(options, 'maxit') % max iteration count: not in algorithm statement
                             % but needed in practice
    if ~isposint(options.maxit)
        error('setdefaults_gradsamp: input "options.maxit" must be a positive integer')
    end
else
    options.maxit = 10000;
end
%
% the following were used in older versions of gradsamp but are no longer valid
%
if isfield(options,'samprad')
    error('setdefaults_gradsamp: options.samprad is no longer valid in GRADSAMP. Use options.samprad0')
end
if isfield(options,'nstart')
    error('setdefaults_gradsamp: options.nstart is no longer valid in GRADSAMP. Only one starting point is allowed')
end
if isfield(options,'cpumax')
    error('setdefaults_gradsamp: options.cpumax is no longer valid in GRADSAMP')
end
if isfield(options,'fvalquit')
    error('setdefaults_gradsamp: options.fvalquit is no longer valid in GRADSAMP')
end
%
% just to avoid any confusion, check that user does not use options.maxit_gradsamp
%
if isfield(options,'maxit_gradsamp')
    error('setdefaults_gradsamp: options.maxit_gradsamp is not valid in GRADSAMP, only in HANSO')
end