function options = setdefaults_bfgs(dimension, options)
%  set check options provided and set defaults for BFGS 
%  to be called only by bfgs

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

% Version note: the main change from version 2.02 to version 2.1 was that
% in version 2.02, the options.nvec default was 0 ONLY if pars.nvar <= 100,
% otherwise it was 10. However, this can be very confusing for a user
% because, especially for nonsmooth problems, limited memory BFGS does not
% perform nearly as well as full BFGS. So, in version 2.1 and later, the
% default is always options.nvec = 0, implying full BFGS. Instead, a
% warning message is printed if the number of variables exceeds 500 and
% options.prtlevel>0.

%
% The following options were moved here from setdefaults.m in the process
% of upgrading to HANSO 3.0 BFGS no longer calls setdefaults.m
%
% if ~isfield(pars, 'nvar')
%    error('setdefaults_bfgs : input "pars" must have a field "nvar" (number of variables)')
% elseif ~isposint(pars.nvar)
%    error('setdefaults_bfgs : input "pars.nvar" (number of variables) must be a positive integer')
% end
% if ~isfield(pars, 'fgname')
%    error('setdefaults_bfgs : input "pars" must have a field "fgname" (name of m-file computing function and gradient)')
% end
if isfield(options, 'maxit')
    % [Modified in 2025:] allow maxit = 0 to test "pure" Gradient Sampling
    if ~isnonnegint(options.maxit)
        error('setdefaults_bfgs : input "options.maxit" must be a nonnegative integer')
    end
else
    options.maxit = 1000;
end
if isfield(options, 'normtol')
    % [Modified in 2025:] allow normtol = 0 since
        % a new stopping condition on obj_tol has been included
    if ~isnonnegreal(options.normtol)
        error('setdefaults_bfgs : input "options.normtol" must be a nonnegative real scalar')
    end
else
    options.normtol = 1.0e-6;
end
if isfield(options, 'fvalquit')
    if ~isreal(options.fvalquit)||~isscalar(options.fvalquit)
        error('setdefaults_bfgs : input "options.fvalquit" must be a real scalar')
    end
else
    options.fvalquit = -inf;
end
if isfield(options, 'xnormquit')
    if ~isreal(options.xnormquit)|~isscalar(options.xnormquit)
        error('setdefaults_bfgs : input "options.xnormquit" must be a real scalar')
    end
else
    options.xnormquit = inf;
end
if isfield(options, 'cpumax')
    if ~isposreal(options.cpumax)
        error('setdefaults_bfgs : input "options.cpumax" must be a positive real scalar')
    end
else
    options.cpumax = inf;
end
if isfield(options, 'prtlevel')
    if options.prtlevel ~= 0 && options.prtlevel ~= 1 && options.prtlevel ~= 2 &&...
            options.prtlevel ~= 3
        error('setdefaults_bfgs : input "options.prtlevel" must be 0,1,2, or 3')
    end
else
    options.prtlevel = 1;
end

if isfield(options, 'print_frequency')
    if ~isposint(options.print_frequency)
        error('setdefaults_bfgs : input "options.print_frequency" must be positive integer')
    end
else
    options.print_frequency = 10000;
end

% line search options
if isfield(options, 'strongwolfe')
    if options.strongwolfe ~= 0 && options.strongwolfe ~= 1
        error('setdefaults_bfgs: input "options.strongwolfe" must be 0 or 1')
    end
else
    % strong Wolfe is very complicated and is bad for nonsmooth functions
    options.strongwolfe = 0;  
end
if isfield(options, 'wolfe1') % Armijo parameter
    % conventionally anything in (0,1), but 0 is usually OK while close to 1 is not
    if ~isreal(options.wolfe1) || options.wolfe1 < 0 || options.wolfe1 > 0.5
        error('setdefaults_bfgs: input "options.wolfe1" must be between 0 and 0.5')
    end
else
    options.wolfe1 = 1e-4; % conventionally this should be positive, and although
                           % zero is usually fine in practice, there are exceptions
end
if isfield(options, 'wolfe2') % Wolfe parameter
    % conventionally should be > wolfe1, but both 0 are OK for e.g. Shor
    if ~isreal(options.wolfe2) || options.wolfe2 < options.wolfe1  || options.wolfe2 >= 1
        error('setdefaults_bfgs: input "options.wolfe2" must be between max(0,options.wolfe1) and 1')
    end
else
    options.wolfe2 = 0.5;  % 0 and 1 are both bad choices
end
if options.strongwolfe
    if options.prtlevel > 0
        if options.wolfe2 > 0
            fprintf('setdefaults_bfgs: strong Wolfe line search selected, but weak Wolfe is usually preferable\n')
            fprintf('(especially if f is nonsmooth)\n')
        else
            fprintf('setdefaults_bfgs: simulating exact line search\n')
        end
    end
    if ~exist('linesch_sw')
        error('"linesch_sw" is not in path: it can be obtained from the NLCG distribution at https://cs.nyu.edu/~overton/software/nlcg/')
    end
else
    if ~exist('linesch_ww')
        error('"linesch_ww" is not in path: it is required for weak Wolfe line search')
    end
end
if isfield(options, 'quitLSfail')
    if options.quitLSfail ~= 0 && options.quitLSfail ~= 1
        error('setdefaults_bfgs: input "options.quitLSfail" must be 0 or 1')
    end
else
    if options.strongwolfe == 1 && options.wolfe2 == 0
        % simulated exact line search, so don't quit if it fails
        options.quitLSfail = 0;
    else  % quit if line search fails
        options.quitLSfail = 1;
    end
end
% other default options
n = dimension;
if isfield(options, 'nvec')
    if ~isnonnegint(options.nvec)
        error('setdefaults_bfgs: input "options.nvec" must be a nonnegative integer')
    end
%%%% elseif n <= 100 THIS WAS IN VERSION 2.02 AND EARLIER BUT IT CAN BE
%%%% VERY CONFUSING FOR A USER: FULL BFGS IS MUCH BETTER WHEN IT IS
%%%% FEASIBLE, ESPECIALLY FOR NONSMOOTH PROBLEMS
else
    options.nvec = 0;  % full BFGS
%%%% else
%%%% options.nvec = 10; % limited memory BFGS
end
if options.nvec == 0 && n > 500 && options.prtlevel > 0
    fprintf('bfgs: when the number of variables is large, consider using the limited memory variant\n');
    fprintf('bfgs: although limited memory may not work as well as full bfgs, especially on nonsmooth problems\n')
    fprintf('bfgs: for limited memory, set options.nvec > 0, for example 10\n');
end 
if isfield(options,'H0')
    % H0 should be positive definite but too expensive to check
    if any(size(options.H0) ~= [n n])
        error('bfgs: input options.H0 must be matrix with order dimension')
    end
    if options.nvec > 0 && ~issparse(options.H0)
        error('bfgs: input "options.H0" must be a sparse matrix when "options.nvec" is positive')
    end
else
    if options.nvec == 0
        options.H0 = eye(n); % identity for full BFGS
    else
        options.H0 = speye(n); % sparse identity for limited memory BFGS
    end
end
if isfield(options, 'scale') % Barzilai-Borwein scaling: see Nocedal&Wright
    if options.scale ~= 0 && options.scale ~= 1
        error('setdefaults_bfgs: input "options.scale" must be 0 or 1')
    end
else
    options.scale = 1; % not important for full BFGS as the initial Hessian
    % is scaled only once, but usually very important for limited memory BFGS
    % where the scaling is done every iteration
end
    
% next option is number of gradients used in termination test
% note: if f is smooth, superlinear convergence will ensure that termination
% takes place before too many gradients are used in the QP optimality check
% so the optimality check will not be expensive in the smooth case
if isfield(options,'ngrad')
    if ~isnonnegint(options.ngrad)
        error('setdefaults_bfgs: input "options.ngrad" must be a nonnegative integer')
    end
else % note this could be more than options.nvec
     % rationale: it is only towards the end that we start accumulating
     % many gradients, and then they may be needed to veryify optimality
    options.ngrad = min([100, 2*dimension, dimension + 10]);
end
if isfield(options,'evaldist') 
    if ~isposreal(options.evaldist) 
        error('setdefaults_bfgs: input "options.evaldist" must be a positive real scalar')
    end
else
    options.evaldist = 1e-4; 
end