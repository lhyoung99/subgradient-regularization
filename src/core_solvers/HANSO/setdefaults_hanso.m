function options = setdefaults_hanso(options)
% check options and set defaults for HANSO version 3.0
% to be called only by hanso
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
% set options
%
% if ~isfield(pars, 'nvar')
%    error('setdefaults_hanso : input "pars" must have a field "nvar" (number of variables)')
% elseif ~isposint(pars.nvar)
%    error('setdefaults_hanso : input "pars.nvar" (number of variables) must be a positive integer')
% end
% if ~isfield(pars, 'fgname')
%    error('setdefaults_hanso : input "pars" must have a field "fgname" (name of m-file computing function and gradient)')
% end
%
% prtlevel applies to BFGS and gradient sampling
%
if isfield(options, 'prtlevel')
    if options.prtlevel ~= 0 && options.prtlevel ~= 1 && options.prtlevel ~= 2 &&...
            options.prtlevel ~= 3
        error('setdefaults_hanso : input "options.prtlevel" must be 0,1,2, or 3')
    end
else
    options.prtlevel = 1;
end

if isfield(options, 'print_frequency')
    if ~isposint(options.print_frequency)
        error('setdefaults_hanso : input "options.print_frequency" must be positive integer')
    end
else
    options.print_frequency = 10000;
end

%
% normtol and evaldist are used for termination in 
% both BFGS and gradient sampling
%
if isfield(options, 'normtol')
    % [Modified in 2025:] allow normtol = 0 since
    % a new stopping condition on obj_tol has been included
    if ~isnonnegreal(options.normtol)
        error('setdefaults_hanso : input "options.normtol" must be a nonnegative scalar')
    end
else
    options.normtol = 1e-4;
end
if isfield(options, 'evaldist')
    if ~isposreal(options.evaldist)
        error('setdefaults_hanso : input "options.evaldist" must be a positive scalar')
    end
else
    options.evaldist = 1e-4;
end

if isfield(options, 'obj_tol')
    if ~isnonnegreal(options.obj_tol)
        error('setdefaults_hanso : input "options.obj_tol" must be a nonnegative scalar')
    end
else
    options.obj_tol = 0;
end
% 
% x0 sets the starting points for (possibly multiple) BFGS runs
%
% if isfield(options,'x0')
%     if size(options.x0,1) ~= pars.nvar
%         error('setdefaults_hanso: input "options.x0" must have "pars.nvar" rows')
%     end
% else
%     options.x0 = randn(pars.nvar, 10); % default is 10 starting points
% end
%
% maxit applies to each BFGS iteration (from each starting point)
%  (maxit_bfgs would be a better name, but keeping maxit for backwards
%  compatibility with HANSO 2.2)
%

% [Modified in 2025:] allow maxit = 0 to test "pure" Gradient Sampling
if isfield(options, 'maxit')
    if ~isnonnegint(options.maxit)
        error('setdefaults_hanso : input "options.maxit" must be a nonnegative integer')
    end
else
    options.maxit = 1000;
end

% [Added in 2025:]
if isfield(options, 'max_oracle_call')
    if ~isnonnegint(options.max_oracle_call)
        error('setdefaults_hanso : input "options.max_oracle_call" must be a nonnegative integer')
    end
else
    options.max_oracle_call = 1000000;
end

% [Modified in 2025:] add time limits
if isfield(options, 'max_time_seconds')
    if ~isposint(options.max_time_seconds)
        error('setdefaults_hanso : input "options.max_time_seconds" must be a positive integer')
    end
else
    options.maxit = 100;
end
%
% nvec > 0 specifies using limited memory BFGS instead of full BFGS
% other BFGS options are set by a call to setdefaultsbfgs from bfgs
%
if isfield(options, 'nvec')
    if ~isnonnegint(options.nvec)
        error('setdefaults_hanso: input "options.nvec" must be a nonnegative integer')
    end
else
    options.nvec = 0;  % full BFGS is much better if number of variables not too large
end
%
% the next three options are only for gradient sampling, not BFGS
% other gradient sampling options are set later by a call to setdefaultGS
% from gradsamp
%
if isfield(options,'maxit_gradsamp')
    if ~isnonnegint(options.maxit_gradsamp)
        error('setdefaults_hanso : input "options.maxit_gradsamp" must be a nonnegative integer')
    end
else
    options.maxit_gradsamp = 0; % default is no gradient sampling
end
if ~isfield(options,'samprad0')
    options.samprad0 = 0.1; % initial sampling radius for gradient sampling
elseif ~isposreal(options.samprad0)
    error('setdefaults_hanso : input "options.samprad0" must be a positive scalar')
end
if ~isfield(options,'target0')
    options.target0 = 0.1; % initial stationarity target for gradient sampling
elseif ~isposreal(options.target0)
    error('setdefaults_hanso : input "options.target0" must be a positive scalar')
end
%
% the following were used in HANSO 2.0, 2.1, 2.2 but are no longer valid
%
if isfield(options,'samprad')
    error('setdefaults_hanso : options.samprad is not valid in HANSO 3.0. Use options.samprad0 instead')
end
if isfield(options,'cpumax')
    error('setdefaults_hanso : options.cpumax is not valid in HANSO 3.0')
end
if isfield(options,'fvalquit')
    error('setdefaults_hanso : options.fvalquit is not valid in HANSO 3.0')
end
%
% the following were used in HANSO 1.0, 1.01 but are no longer valid
%
if isfield(options, 'phasemaxit')
    error('setdefaults_hanso : options.phasemaxit is not valid in HANSO 3.0')
end
if isfield(options, 'phasenum')
    error('setdefaults_hanso : options.phasenum is not valid in HANSO 3.0')
end