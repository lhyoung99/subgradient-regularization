function results = bfgs1run(x0, func_value, subgrad, flags, options)
%{
Make a single run of BFGS from one starting point
To be called only by bfgs.m

Outputs: 
  A struct containing:
   x: final iterate
   f: final function value
   d: final smallest vector in convex hull of saved gradients
   H: final inverse Hessian approximation
   iter: number of iterations
   info: reason for termination
    0: tolerance on smallest vector in convex hull of saved gradients met
    1: max number of iterations reached
    2: f reached target value
    3: norm(x) exceeded limit
    4: cpu time exceeded limit
    5: f or g is inf or nan at initial point
    6: direction not a descent direction (because of rounding)
    7: line search bracketed minimizer but Wolfe conditions not satisfied
    8: line search did not bracket minimizer: f may be unbounded below 
   X: iterates where saved gradients were evaluated (with X(:,1) = x)
   G: gradients evaluated at these points
   w: weights defining convex combination d = G*w
==== [Dropped in 2025/02/17] ====
   fevalrec: record of all function evaluations in the line searches
   xrec: record of x iterates
   Hrec: record of H iterates
=================================
  Send comments/bug reports to Michael Overton, mo1@nyu.edu.
  For HANSO 3.0, 2021 (with only cosmetic changes to version 2.02, 2013)
%}
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

% Initial parameters
n = length(x0); 
normtol = options.normtol;
obj_tol = options.obj_tol;
fvalquit = options.fvalquit;
xnormquit = options.xnormquit;
% cpufinish = cputime + options.max_time_seconds;
maxit = options.maxit;
nvec = options.nvec;
prtlevel = options.prtlevel;
print_frequency = options.print_frequency;
strongwolfe = options.strongwolfe;
wolfe1 = options.wolfe1;
wolfe2 = options.wolfe2;
quitLSfail = options.quitLSfail;
ngrad = options.ngrad;
evaldist = options.evaldist;

% Extract opt_val & opt_soln from 'flags' if exist
opt_val = [];
if isfield(flags, 'opt_val')
    opt_val = flags.opt_val;
end

opt_soln = [];
if isfield(flags, 'opt_soln')
    opt_soln = flags.opt_soln;
end

% Initialize iterate
H0 = options.H0;
H = H0; % sparse for limited memory BFGS 
scale = options.scale;
x = x0;
f = func_value(x);
g = subgrad(x);
d = g;
G = g;
X = x; 
nG = 1;
w = 1;
dnorm = norm(g);
oracle = 1;

% % so outputs defined if quit immediately
% fevalrec{1} = nan; % cell array
% xrec = nan*ones(n,1); % not cell array
% Hrec{1} = nan; % cell array

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
if nvec > 0 % limited memory BFGS
    S = [];
    Y = [];
end
iter = 0;

% Record data 
elapsed_time_list(end+1) = toc(start_time);
num_oracle_list(end+1) = oracle;
history = [history, x(:)];
distance_list(end+1) = distance_or_empty(distance_current);
objective_values_list(end+1) = f;

if isnaninf(f) % better not to generate an error return
    if prtlevel > 0
        fprintf('bfgs: f is infinite or nan at initial iterate\n')
    end
    info = 5;
    results = false; 
    return
elseif isnaninf(g)
    if prtlevel > 0
        fprintf('bfgs: gradient is infinite or nan at initial iterate\n')
    end
    info = 5;
    results = false;
    return
elseif dnorm < normtol
    if prtlevel > 0
        fprintf('bfgs: tolerance on gradient satisfied at initial iterate\n')
    end
    info = 0;
    results = compile_statistics(...
            x, f, d, H, iter, info, X, G, w, ...%fevalrec, xrec, Hrec, ...
            num_oracle_list, elapsed_time_list, objective_values_list, distance_list, history ...
            );
    return
elseif ~isempty(opt_val) && f - opt_val <= obj_tol
    if prtlevel > 0
        fprintf('bfgs: tolerance on obj val satisfied at initial iterate\n')
    end
    info = 2;
    results = compile_statistics(...
            x, f, d, H, iter, info, X, G, w, ...%fevalrec, xrec, Hrec, ...
            num_oracle_list, elapsed_time_list, objective_values_list, distance_list, history ...
            );
    return
elseif f < fvalquit
    if prtlevel > 0
        fprintf('bfgs: below target objective at initial iterate\n')
    end
    info = 2;
    results = compile_statistics(...
            x, f, d, H, iter, info, X, G, w, ...%fevalrec, xrec, Hrec, ...
            num_oracle_list, elapsed_time_list, objective_values_list, distance_list, history ...
            );
    return
elseif norm(x) > xnormquit
    if prtlevel > 0
        fprintf('bfgs: norm(x) exceeds specified limit at initial iterate\n')
    end
    info = 3;
    results = compile_statistics(...
            x, f, d, H, iter, info, X, G, w, ...%fevalrec, xrec, Hrec, ...
            num_oracle_list, elapsed_time_list, objective_values_list, distance_list, history ...
            );
    return
end

for iter = 1:maxit
    if nvec == 0 % full BFGS
        p = -H*g;
    else % limited memory BFGS
        p = -hgprod(H, g, S, Y);  % not H0, as in previous version
    end
    gtp = g'*p;
    if gtp >= 0 || isnan(gtp) % in rare cases, H could contain nans
       if prtlevel > 0
          fprintf('bfgs: not descent direction, quit at iteration %d, f = %10.7e, dnorm = %5.1e\n',...
              iter, f, dnorm)
       end
       info = 6;
       results = compile_statistics(...
                x, f, d, H, iter, info, X, G, w, ...%fevalrec, xrec, Hrec, ...
                num_oracle_list, elapsed_time_list, objective_values_list, distance_list, history ...
                );
       return
    end
    gprev = g;  % for BFGS update
    if ~strongwolfe % weak Wolfe line search is the default
        [alpha, x, f, g, fail, ~, ~, fevalrecline, ls_oracle] = ...
                linesch_ww(x, f, g, p, func_value, subgrad, wolfe1, wolfe2, fvalquit, prtlevel);
        oracle = oracle + ls_oracle;
        % Update distance to opt_soln
        if ~isempty(opt_soln)
            distance_current = norm(x - opt_soln);
        end

        % Record data
        elapsed_time_list(end+1) = toc(start_time);
        num_oracle_list(end+1) = oracle;
        history = [history, x(:)];
        distance_list(end+1) = distance_or_empty(distance_current);
        objective_values_list(end+1) = f;
    end
    % for the optimality check:
    % discard the saved gradients iff the new point x is not sufficiently
    % close to the previous point and replace them by new gradient 
    if alpha*norm(p) > evaldist
        nG = 1;
        G = g;
        X = x;
    % otherwise add new gradient to set of saved gradients, 
    % discarding oldest if already have ngrad saved gradients
    elseif nG < ngrad
        nG = nG + 1;
        G =  [g G];
        X = [x X];
    else % nG = ngrad
        G = [g G(:,1:ngrad-1)];
        X = [x X(:,1:ngrad-1)];
    end
    % optimality check: compute smallest vector in convex hull of qualifying 
    % gradients: reduces to norm of latest gradient if ngrad == 1, and the
    % set must always have at least one gradient: could gain efficiency
    % here by updating previous QP solution
    if nG > 1
        [w,d] = qpspecial(G); % Anders Skajaa code for this special QP
    else
        w = 1;
        d = g; 
    end
    dnorm = norm(d);

    % xrec(:,iter) = x;
    % fevalrec{iter} = fevalrecline; % function vals computed in line search
    % Hrec{iter} = H;

    if prtlevel > 1 && mod(iter, print_frequency) == 0
        nfeval = length(fevalrecline);
        fprintf('bfgs: iter %d: nfevals = %d, step = %5.1e, f = %10.7e, nG = %d, dnorm = %5.1e\n', ...
            iter, nfeval, alpha, f, nG, dnorm)
    end
    if f < fvalquit % this is checked inside the line search
        if prtlevel > 0
            fprintf('bfgs: reached target objective, quit at iteration %d \n', iter)
        end
        info = 2;
        results = compile_statistics(...
                x, f, d, H, iter, info, X, G, w, ...%fevalrec, xrec, Hrec, ...
                num_oracle_list, elapsed_time_list, objective_values_list, distance_list, history ...
                );
        return
    elseif norm(x) > xnormquit % this is not checked inside the line search
        if prtlevel > 0
            fprintf('bfgs: norm(x) exceeds specified limit, quit at iteration %d \n', iter)
        end
        info = 3;
        results = compile_statistics(...
                x, f, d, H, iter, info, X, G, w, ...%fevalrec, xrec, Hrec, ...
                num_oracle_list, elapsed_time_list, objective_values_list, distance_list, history ...
                );
        return
    end
    if fail == 1 % line search failed (Wolfe conditions not both satisfied)
        if ~quitLSfail
            if prtlevel > 1
                fprintf('bfgs: continue although line search failed\n')
            end
        else % quit since line search failed
            if prtlevel > 0
                fprintf('bfgs: quit at iteration %d, f = %10.7e, dnorm = %5.1e\n', iter, f, dnorm)
            end
            info = 7;
            results = compile_statistics(...
                    x, f, d, H, iter, info, X, G, w, ...%fevalrec, xrec, Hrec, ...
                    num_oracle_list, elapsed_time_list, objective_values_list, distance_list, history ...
                    );
            return
        end
    elseif fail == -1 % function apparently unbounded below
        if prtlevel > 0
           fprintf('bfgs: f may be unbounded below, quit at iteration %d, f = %10.7e\n', iter, f)
        end
        info = 8;
        results = compile_statistics(...
                x, f, d, H, iter, info, X, G, w, ...%fevalrec, xrec, Hrec, ...
                num_oracle_list, elapsed_time_list, objective_values_list, distance_list, history ...
                );
        return
    end

    %% Termination condition
    if dnorm <= normtol || (~isempty(opt_val) && f - opt_val <= obj_tol)
        if prtlevel > 0 
            if dnorm <= normtol && nG == 1 
                fprintf('bfgs: gradient norm below tolerance, quit at iteration %d, f = %10.7e\n', iter, f')
            elseif dnorm <= normtol
                fprintf('bfgs: norm of smallest vector in convex hull of gradients below tolerance, quit at iteration %d, f = %10.7e\n', iter, f')
            else
                fprintf('bfgs: obj val below tolerance, quit at iteration %d, f = %10.7e\n', iter, f')
            end
        end
        info = 0;
        results = compile_statistics(...
                x, f, d, H, iter, info, X, G, w, ...%fevalrec, xrec, Hrec, ...
                num_oracle_list, elapsed_time_list, objective_values_list, distance_list, history ...
                );
        return
    end

    current_time = toc(start_time);
    if current_time > options.max_time_seconds
        if prtlevel > 0
            fprintf('bfgs: cpu time limit exceeded, quit at iteration %d\n', iter)
        end
        info = 4;
        results = compile_statistics(...
                x, f, d, H, iter, info, X, G, w, ...%fevalrec, xrec, Hrec, ...
                num_oracle_list, elapsed_time_list, objective_values_list, distance_list, history ...
                );
        return
    end

    if oracle > options.max_oracle_call
        info = 4;
        if prtlevel > 0
            fprintf('bfgs: max number of oracle calls exceeded\n')
        end
        results = compile_statistics(...
                x, f, d, H, iter, info, X, G, w, ...%fevalrec, xrec, Hrec, ...
                num_oracle_list, elapsed_time_list, objective_values_list, distance_list, history ...
                );
        return
    end

    s = alpha*p;
    y = g - gprev;
    sty = s'*y;    % successful line search ensures this is positive
    if nvec == 0   % perform rank two BFGS update to the inverse Hessian H
        if sty > 0 
            if iter == 1 && scale
                % for full BFGS, Nocedal and Wright recommend scaling before 
                % the first update only
                H = (sty/(y'*y))*H; 
            end
            % for formula, see Nocedal and Wright's book
            %M = I - rho*s*y', H = M*H*M' + rho*s*s', so we have
            %H = H - rho*s*y'*H - rho*H*y*s' + rho^2*s*y'*H*y*s' + rho*s*s'
            % note that the last two terms combine: (rho^2*y'Hy + rho)ss'
            rho = 1/sty;
            Hy = H*y;
            rhoHyst = (rho*Hy)*s';  
            % old version: update may not be symmetric because of rounding 
            % H = H - rhoHyst' - rhoHyst + rho*s*(y'*rhoHyst) + rho*s*s';  
            % new in version 2.02: make H explicitly symmetric
            % also saves one outer product
            % in practice, makes little difference, except H=H' exactly
            ytHy = y'*Hy; % could be < 0 if H not numerically pos def
            sstfactor = max([rho*rho*ytHy + rho,  0]);
            sscaled = sqrt(sstfactor)*s;
            H = H - (rhoHyst' + rhoHyst) + sscaled*sscaled';
            % alternatively add the update terms together first: does
            % not seem to make significant difference
            % update = sscaled*sscaled' - (rhoHyst' + rhoHyst);
            % H = H + update;
        else % should not happen unless line search fails, and in that case should normally have quit
            if prtlevel > 1
                fprintf('bfgs: sty <= 0, skipping BFGS update at iteration %d \n', iter)
            end
        end
    else % save s and y vectors for limited memory update
        s = alpha*p;
        y = g - gprev;
        if iter <= nvec
            S = [S s];
            Y = [Y y];
        else % could be more efficient here by avoiding moving the columns
            S = [S(:,2:nvec) s];
            Y = [Y(:,2:nvec) y];
        end
        if scale 
            H = ((s'*y)/(y'*y))*H0;  % recommended by Nocedal-Wright
        end
    end
end % for loop
iter = maxit;
if prtlevel > 0
    fprintf('bfgs: %d iterations reached, f = %10.7e, dnorm = %5.1e\n', maxit, f, dnorm)
end
info = 1; % quit since max iterations reached
results = compile_statistics(...
        x, f, d, H, iter, info, X, G, w, ...%fevalrec, xrec, Hrec, ...
        num_oracle_list, elapsed_time_list, objective_values_list, distance_list, history ...
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
    x, f, d, H, iter, info, X, G, w, ... % fevalrec, xrec, Hrec, ...
    num_oracle_list, elapsed_time_list, objective_values_list, ...
    distance_list, history ...
)
    results_struct.x = x;
    results_struct.f = f;
    results_struct.d = d;
    results_struct.H = (H + H')/2;
    results_struct.iter = iter;
    results_struct.info = info;
    results_struct.X = X;
    results_struct.G = G;
    results_struct.w = w;
    % results_struct.fevalrec = fevalrec;
    % results_struct.xrec = xrec;
    % results_struct.Hrec = Hrec;
    results_struct.num_oracle_list = num_oracle_list;
    results_struct.elapsed_time_list = elapsed_time_list;
    results_struct.objective_values_list = objective_values_list;
    results_struct.distance_list = distance_list;
    results_struct.history = history;
end
