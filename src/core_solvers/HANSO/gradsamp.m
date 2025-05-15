function results = gradsamp(x0, func_value, subgrad, flags, options)
%GRADSAMP: Gradient Sampling, 2021 
%  Gradient sampling algorithm for nonsmooth, nonconvex minimization.
%  Intended for nonconvex functions that are continuous everywhere and for 
%  which the gradient can be computed at most points, but which are known 
%  to be nondifferentiable at some points, typically including minimizers.
%  Completely rewritten in 2021 to implement Algorithm GS in our survey paper:
%  [1] J.V. Burke, F.E. Curtis, A.S. Lewis, M.L. Overton and L.E.A. Simões, 
%    Gradient Sampling Methods for Nonsmooth Optimization,
%    In: Numerical Nonsmooth Optimization, edited by A. Bagirov et al, 
%    Springer (2020), pp. 201-225, https://arxiv.org/abs/1804.11003.
%  See below for more details. See [1] for convergence properties.
%
%   Calls:  [x, f, d, X, G, w, termcode] = gradsamp(func_value,subgrad, options) 
%    or version with more output args, see function header
%
%   Input parameters
%    func_value
%    subgrad
%    options is an optional struct 
%       options.eps0: initial value of sampling radius 
%          (default: 0.1)
%       options.nu0: initial stationarity target
%          (default: 0.1)
%       options.ngrad: number of sampled gradients per iterate
%          normally should be >= nvar+1, but can be set smaller,
%          even to 0 to implement steepest descent
%          (default: 2*nvar) 
%       options.beta: Armijo-like parameter used in backtracking line search
%          (default: 1e-4)
%       options.gamma: contraction factor in backtracking line search
%          (default: 0.5)
%       options.eps_opt: sampling radius termination tolerance
%          (default: 1e-6)
%       options.nu_opt: stationarity termination tolerance
%          (default: 1e-6)
%       options.theta_eps: reduction factor for sampling radius
%          (default: 0.1)
%       options.theta_nu: reduction factor for stationarity target
%          (default: 0.1)
%       options.delta: use to terminate line search if t gets too small
%          (default: 1e-16)
%       options.maxit: max number of total iterations (not per sampling radius)
%          (default: 10000)
%       options.prtlevel: print level, between 0 (none) and 3 (verbose).
%          (default: 1)
%      Note: the last 3 options are implementation details and do not
%        appear in the definition of Algorithm GS in [1].
% ===============[ New options added in 2025 ]================
%       options.print_frequency: the gap of printing k
%          (default: 10000)
%       options.obj_tol: tolerance of objective value for stopping
%          (default: empty)
% ============================================================
%   Output parameters 
%    A struct containing:
%    x: final value of iterate x
%    f: function value at x
%    d: the smallest vector in the convex hull of the gradient at x and 
%       the sampled gradients near x [somewhat confusingly, this is called 
%       "g" in Algorithm GS in [1], but it is not a gradient]
%       The norm of d is a measure of approximate Clarke stationarity
%    X: X(:,1) is the same as x, and X(:,j), j=2,3,...,options.ngrad+1 are 
%       points where the gradients were last evaluated (sampled)
%    G: G(:,j) is the gradient of the function at X(:,j) 
%    w: w is the vector of nonnegative weights summing to one such that
%       d = G*w  (length options.ngrad + 1)
%    termcode: termination code:
%       0: satisfied stopping criterion
%       1: -d is not a descent direction
%       2: line search failure
%       3: maximum iterations/time reached
%    frec: record of function values per iteration
%    dnormrec: record of ||d|| per iteration
%    epsrec: record of sampling radii per iteration
%    evalsrec: record of number of function values computed in line search, 
%              per iteration
%    xrec: record of all iterates (except those computed in line search)
%    flinerec: cell array of all computed function values including those 
%              in line search
%  Note: the differentiability check in lines 8-10 of Algorithm GS is
%  omitted, since there is really no way to implement it in floating point.
%  Otherwise this implementation should be exactly as described in the
%  survey paper [1] (with some additional options parameters as noted above).
%
%  This version of GRADSAMP is incorporated in HANSO version 3.0, but it
%  can be called as a stand-alone function. 

%  As explained in [1], Algorithm GS is based on a 2007 paper of Kiwiel
%  which suggested several modifications to the original algorithm in: 
%  [2] J.V. Burke, A.S. Lewis and M.L. Overton, A Robust Gradient Sampling 
%  Algorithm for Nonsmooth, Nonconvex Optimization
%  SIAM J. Optimization 15 (2005), pp. 751-779.
%  The older versions of GRADSAMP in HANSO versions 1.0-2.2 were based 
%  on [2] but they included a number of changes to the code for which 
%  results were reported in [2], in particular:
%  (a) reverting to the unnormalized search direction we used in 2002 as 
%  later also recommended by Kiwiel [retained in this version], 
%  (b) using our standard "weak Wolfe" line search code linesch_ww, with the 
%  Armijo and Wolfe parameters both set to zero, instead of the simpler line 
%  search in the BLO paper [no longer used in this version], and 
%  (c) using Andres Skajaa's "qpspecial" QP solver instead of quadprog to 
%  compute the search direction [retained in this version, with tightened 
%  tolerances, after some testing against quadprog and Curtis' implementation 
%  of Kiwiel's QP algorithm, both of which are intended for more general QPs].
%
%  For the original code used in the experiments reported in [2], 
%  see www.cs.nyu.edu/overton/papers/gradsamp/alg/
%
%  Send comments/bug reports to Michael Overton, mo1@nyu.edu.

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
%% [Dropped in 2025/02/17 by wrapping with a construct]
% if nargout > 8
%     getrecs = true; % get record vectors
%     % these records are values at the end of the for loop, but it might
%     % terminate in the middle of the loop in which case these values will
%     % not be updated at the end of the loop, and could be empty if
%     % termination takes place on the first iteration
%     frec = []; 
%     dnormrec = []; 
%     epsrec = [];
%     evalsrec = [];
%     xrec = [];
%     flinerec = {};
% else
%     getrecs = false;
% end
%
% check pars and options fields and set defaults for options
%

% Initial parameters
options = setdefaults_gradsamp(flags.nb_parameters, options); 
ngrad = options.ngrad; % number of gradients to sample
eps = options.eps0; % sampling radius
nu = options.nu0; % stationarity target
eps_opt = options.eps_opt; % sampling radius termination tolerance
nu_opt = options.nu_opt; % stationarity termination tolerance
theta_eps = options.theta_eps; % reduction factor for sampling radius
theta_nu = options.theta_nu; % reduction factor for stationarity target
maxit = options.maxit; % max number of iterations
prtlevel = options.prtlevel; % print level
print_frequency = options.print_frequency;
obj_tol = options.obj_tol;

opt_soln = [];
if isfield(flags, 'opt_soln')
    opt_soln = flags.opt_soln;
end
% cpufinish = cputime + options.max_time_seconds;

% Initialize iterate
f0 = func_value(x0);
g0 = subgrad(x0);
f = f0;
x = x0;
g = g0;
num_oracle = 1;

num_oracle_list = [];
objective_values_list = [];
distance_list = [];
elapsed_time_list = [];
history = [];
if ~isempty(opt_soln)
    distance_current = norm(x - opt_soln);
else
    distance_current = [];
end

% Timer start
start_time = tic;
if isnan(f0)
    error('gradsamp: function is NaN at initial point\n')
elseif f0 == inf 
    error('gradsamp: function is infinite at initial point\n')
end

% Record data
elapsed_time_list(end+1) = toc(start_time);
num_oracle_list(end+1) = num_oracle;
history = [history, x(:)];
distance_list(end+1) = distance_or_empty(distance_current);
objective_values_list(end+1) = f;

%
% there is only one loop in Algorithm GS, not counting the line search
%
% Line 1 of Algorithm GS (typo: should be "for k in 0,1,2,...", not k in N)
%
for k = 0:maxit-1 % maxit is the max number of iterations, including k=0
    %
    % Line 2 of Algorithm GS:
    % sample gradients of f uniformly from ball with radius eps around x
    %
    [X,G,oracle] = samplegrads(x,eps,ngrad,func_value,subgrad);
    num_oracle = num_oracle + oracle;

    % Record data
    elapsed_time_list(end+1) = toc(start_time);
    num_oracle_list(end+1) = num_oracle;
    history = [history, x(:)];
    distance_list(end+1) = distance_or_empty(distance_current);
    objective_values_list(end+1) = f;
    %
    % Line 3 of Algorithm GS:
    % solve QP subproblem: d is called g in Algorithm GS
    % solve with Skajaa's qp_special, modified to have tighter tolerances
    % testing shows this is just as accurate and much faster than quadprog
    %
    % prepend x to X and g to G: in previous versions, this was not done,
    % as the convention was that qpspecial was called with argument [g G]
    %
    X = [x X]; % for output argument
    G = [g G]; % for passing to qpspecial and output argument
    [w,d] = qpspecial(G); % get smallest vector d in convex hull of columns
    dnorm = norm(d);
    %
    %% Line 4 of Algorithm GS: termination condition
    %
    if (dnorm <= nu_opt && eps <= (1+1e-14)*eps_opt) || (~isempty(flags.opt_val) && f - flags.opt_val <= obj_tol)
        % Algorithm GS says dnorm <= nu_opt and eps <= eps_opt, but
        % testing eps <= eps_opt may result in one extra reduction
        % when using powers of 10, not powers of 2, for eps0 and eps_opt,
        % because of rounding errors when e.g. repeatedly multiplying by 0.1
        % (this care is not necessary when comparing dnorm to nu_opt because
        % dnorm is computed from solving a QP, not from a fixed parameter
        % or repeated multiplications of a parameter by a constant factor)
        % 
        termcode = 0; % successful termination
        if prtlevel > 0
            fprintf('gradsamp: k=%d, terminate as ||d|| is %g and eps is %g, f = %10.7e\n', k, dnorm,eps,f)
        end
        results = compile_statistics(...
                x, f, d, X, G, w, termcode, ...%frec, dnormrec, epsrec, evalsrec, xrec, flinerec, ...
                num_oracle_list, elapsed_time_list, objective_values_list, ...
                distance_list, history ...
                );
        return
    end
    
    current_time = toc(start_time);
    if current_time > options.max_time_seconds
        termcode = 3;
        if prtlevel > 0
            fprintf('gradsamp: cpu time limit exceeded\n')
        end
        results = compile_statistics(...
                x, f, d, X, G, w, termcode, ...%frec, dnormrec, epsrec, evalsrec, xrec, flinerec, ...
                num_oracle_list, elapsed_time_list, objective_values_list, ...
                distance_list, history ...
                );
        return
    end
    
    if num_oracle > options.max_oracle_call
        termcode = 3;
        if prtlevel > 0
            fprintf('gradsamp: max number of oracle calls exceeded\n')
        end
        results = compile_statistics(...
                x, f, d, X, G, w, termcode, ...%frec, dnormrec, epsrec, evalsrec, xrec, flinerec, ...
                num_oracle_list, elapsed_time_list, objective_values_list, ...
                distance_list, history ...
                );
        return
    end
    %
    %% Line 5 of Algorithm GS: check if current stationarity target met
    %
    if dnorm <= nu % not nu_opt! 
        % although, it sometimes works fine to use nu_opt here, and indeed
        % in the BLO 2005 paper [2], that is what we did: see p.767, para.4
        %
        %% Line 6 of Algorithm GS: reduce stationarity target and sampling radius
        %
        nu = theta_nu*nu;  
        eps = theta_eps*eps;
        % no change in x (equivalently, t=0 in Algorithm GS)
        xnew = x;
        fnew = f;
        gnew = g;
        evals = 0;
        fline = [];
        if prtlevel > 0 && mod(k, print_frequency) == 0
            fprintf('gradsamp: k=%d, ||d||=%g, reducing nu to %g and eps to %g, f = %10.7e\n',k,dnorm,nu,eps,f)
        end
    else
        % stationarity target not met, so use -d for line search
        % in theory, -d should be a descent direction
        gtd = g'*d;
        if gtd <= 0 %-1e-7 (for more iterations) % g'(-d) >= 0
            % this test is not in Algorithm GS because it cannot happen in
            % exact arithmetic. However, if gtd <= 0 because of rounding 
            % error, -d is not a descent direction, so terminate algorithm.
            % [Might instead consider reducing eps and continuing; 
            % this is what we did in the BLO 2005 paper [2]: see p.767, para.5,
            % but this would undermine the nice theory of Algorithm GS.]
            termcode = 1;
            if prtlevel > 0
                fprintf('gradsamp: k=%d, terminate as -d is not a descent direction, g''d = %g, f=%10.7e\n',k,gtd,f)
            end
            results = compile_statistics(...
                x, f, d, X, G, w, termcode, ...%frec, dnormrec, epsrec, evalsrec, xrec, flinerec, ...
                num_oracle_list, elapsed_time_list, objective_values_list, ...
                distance_list, history ...
                );
            return
        end
        %
        %% Line 7 of Algorithm GS: backtracking line search along 
        % direction -d, leaving nu and eps unchanged
        %
        [xnew,fnew,gnew,evals,ls_termcode,fline] = linesch_gradsamp(x,f,-d,func_value,subgrad,options);
        num_oracle = num_oracle + evals;
        % Update distance to opt_soln
        if ~isempty(opt_soln)
            distance_current = norm(xnew - opt_soln);
        end

        % Record data
        elapsed_time_list(end+1) = toc(start_time);
        num_oracle_list(end+1) = num_oracle;
        history = [history, xnew(:)];
        distance_list(end+1) = distance_or_empty(distance_current);
        objective_values_list(end+1) = fnew;

        if prtlevel > 1 && mod(k, print_frequency) == 0
            fprintf('gradsamp: k=%d, nu=%g, eps=%g,||d||=%12.6e, line search: evals=%d, ls_termcode=%d, fnew=%10.7e\n',...
                k,nu,eps,dnorm,evals,ls_termcode,fnew)
        end
        if ls_termcode > 0 
            % this test is not in Algorithm GS because it cannot happen in 
            % exact arithmetic. However, if the line search fails due to
            % rounding errors, terminate the algorithm.
            % [Might instead consider reducing eps and continuing; 
            % this is what we did in the BLO 2005 paper: see p.767, para.5,
            % but this would undermine the nice theory of Algorithm GS.]
            termcode = 2;
            if prtlevel > 0
                fprintf('gradsamp: k=%d, terminate since line search failed, f=%10.7e\n',k,f)
            end
            results = compile_statistics(...
                x, f, d, X, G, w, termcode, ...%frec, dnormrec, epsrec, evalsrec, xrec, flinerec, ...
                num_oracle_list, elapsed_time_list, objective_values_list, ...
                distance_list, history ...
                );
            return
        end
    end % end of else part
    %
    %% Lines 8-10 of Algorithm GS: update x, but skip the differentiability 
    %  check, which makes no sense in finite precision
    %
    x = xnew; % unchanged if stationarity target and sampling radius were reduced
    f = fnew;
    g = gnew;
    %
    % if requested, report records of each function value (not including
    % those evaluated in the line search), norm of search direction,
    % sampling radius, number of function values computed in the line 
    % search, and the x iterates (again not including those evaluated in
    % the line search). Note that the final values are not included if we 
    % already broke out of the loop earlier, but this is not important 
    % since the important output args f, d and x contain final values.
    %
    % if getrecs 
    %     kp1 = k+1; % since k starts at 0
    %     frec(kp1) = f; 
    %     dnormrec(kp1) = dnorm; 
    %     epsrec(kp1) = eps;
    %     evalsrec(kp1) = evals;
    %     xrec(:,kp1) = x;
    %     flinerec{kp1} = fline; % function vals computed in line search
    % end
end % end of for loop
%
% reached maxit iterations: there is no such possibility in Algorithm GS, 
% which is an infinite loop
%
termcode = 3;
if prtlevel > 0
    fprintf('gradsamp: %d iterations reached, f=%10.7e\n',maxit,f)
end
%
% Repeat Lines 2 and 3 of Alg GS since maxit iterations were reached
% and either x or eps was updated: otherwise, either x will not be the
% first column of X or the gradients will not be sampled with radius
% eps. Effectively this means executing the first part of the loop one
% more time. However, do not update dnormrec since it reports the value
% of ||d|| at the end of each iteration.
%
[X,G,oracle] = samplegrads(x,eps,ngrad,func_value,subgrad);
X = [x X]; 
G = [g G]; 
num_oracle = num_oracle + oracle;
[w,d] = qpspecial(G); 
results = compile_statistics(...
        x, f, d, X, G, w, termcode, ...%frec, dnormrec, epsrec, evalsrec, xrec, flinerec, ...
        num_oracle_list, elapsed_time_list, objective_values_list, ...
        distance_list, history ...
        );
end

%% ===================== Helper functions =====================
function val = distance_or_empty(distance_current)
    if isempty(distance_current)
        val = NaN;
    else
        val = distance_current;
    end
end

function results_struct = compile_statistics(...
    x, f, d, X, G, w, termcode, ... %frec, dnormrec, epsrec, evalsrec, xrec, flinerec, ...
    num_oracle_list, elapsed_time_list, objective_values_list, ...
    distance_list, history ...
)
    results_struct.x = x;
    results_struct.f = f;
    results_struct.d = d;
    results_struct.X = X;
    results_struct.G = G;
    results_struct.w = w;
    results_struct.termcode = termcode;
    % results_struct.frec = frec;
    % results_struct.dnormrec = dnormrec;
    % results_struct.epsrec = epsrec;
    % results_struct.evalsrec = evalsrec;
    % results_struct.xrec = xrec;
    % results_struct.flinerec = flinerec;
    results_struct.num_oracle_list = num_oracle_list;
    results_struct.elapsed_time_list = elapsed_time_list;
    results_struct.objective_values_list = objective_values_list;
    results_struct.distance_list = distance_list;
    results_struct.history = history;
end
