%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% optval_LICQ.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [func_value, subgrad, G_oracle, flags] = optval_LICQ(dimension_x, dimension_y, dimension_y_cons, seed)
% OPTVAL_LICQ
%   f(x) = ||x||^4 + min{ -(c + D * x)'y + 0.5 * y' * Q * y }  subject to A*x + W*y <= b
%
% Outputs:
%   func_value : function handle for f(x)
%   subgrad    : function handle optval_LICQ_subgrad(x) that returns the exact gradient
%   G_oracle   : function handle optval_LICQ_G_oracle(x, epsilon) that approximates the subgradient
%   flags      : structure with additional information
    if nargin < 1, dimension_x = 10; end
    if nargin < 2, dimension_y = 5; end
    if nargin < 3, dimension_y_cons = 5; end
    if nargin < 4, seed = 3407; end
    rng(seed);
    n = dimension_x; m = dimension_y; r = dimension_y_cons;
    A = randn(r, n) / sqrt(r);
    b = randn(r, 1) / sqrt(r);
    c = randn(m, 1) / sqrt(m);
    D = randn(m, n) / sqrt(m);
    tmp = randn(r, m);
    [QO, ~] = qr(tmp');
    W = QO(:,1:r)'; 
    tmp2 = randn(m, m) / sqrt(m);
    Q = tmp2' * tmp2;
    
    func_value = @(x) optval_LICQ_func_value(x, A, b, c, D, W, Q);
    G_oracle   = @(x, epsilon) optval_LICQ_G_oracle(x, epsilon, A, b, c, D, W, Q);
    subgrad    = @(x) optval_LICQ_subgrad(x, A, b, c, D, W, Q);
    
    flags.function_name = 'opt_val_function_LICQ';
    flags.opt_val = [];
    flags.opt_soln = [];
    flags.nb_parameters = n;
    flags.params_for_legend = [];
end

function val = optval_LICQ_func_value(x, A, b, c, D, W, Q)
    m = length(c);
    % %% (Drop due to inefficiency) Call Gurobi to solve QP %%
    % model.A = sparse(W);
    % model.rhs = b - A*x;
    % model.sense = '<';
    % model.lb = -inf(n_y, 1);
    % model.obj = (c + D*x);
    % model.Q = sparse(-Q/2);
    % model.modelsense = 'max';
    % params.OutputFlag = 0;
    % result = gurobi(model, params);
    % if strcmp(result.status, 'OPTIMAL')
    %     val = result.objval + 50 * norm(x)^2;
    %     soln = result.x;
    % else
    %     val = [];
    %     error('Optimal value function does not have a solution for current x');
    % end

    %% Call Matlab Opt Toolbox to solve QP %%
    options.Display = 'off';
    options.Algorithm = 'active-set';  % 'interior-point-convex' (default)/'trust-region-reflective'/'active-set'
    % For a 'trust-region-reflective' equality-constrained problem, the default value is 1e-6.
    % For a 'trust-region-reflective' bound-constrained problem, the default value is 100*eps, about 2.2204e-14.
    % For the 'interior-point-convex' and 'active-set' algorithms, the default value is 1e-8.
    % options.OptimalityTolerance = 1e-14; 
    
    y0 = zeros(m,1);
    lhs = [W; -W];
    rhs = [b - A*x; -b + A*x + ones(m,1)];
    [~, fval, exitflag] = quadprog(Q, -(c + D*x), lhs, rhs, [], [], [], [], y0, options);
    % [~, fval, exitflag] = linprog(-(c + D*x), W, b - A*x, [], [], [], [], options);
    if exitflag == -3
        error('Optimal value function does not have a solution for current x');
    else
        val = - fval + norm(x)^4;
    end
end

function approx_grad = optval_LICQ_G_oracle(x, epsilon, A, b, c, D, W, Q)
    m = length(c);
    % %% (Drop due to inefficiency) Call Gurobi to solve QP %%
    % model1.A = sparse(W);
    % model1.rhs = b - A*x;
    % model1.sense = '<';
    % model1.lb = -inf(m, 1);
    % model1.obj = (c + D*x);
    % model1.Q = sparse(-Q/2);
    % model1.modelsense = 'max';
    % params.OutputFlag = 0;
    % result1 = gurobi(model1, params);
    % if ~strcmp(result1.status, 'OPTIMAL')
    %     error('First QP not solved to optimal.');
    % end
    % duals = result1.pi;
    
    %% Call Matlab Opt Toolbox to solve QP %%
    options.Display = 'off';
    % options.Algorithm = 'active-set';  % 'interior-point-convex' (default)/'trust-region-reflective'/'active-set'
    
    % y0 = zeros(m,1);
    lhs = [W; -W];
    rhs = [b - A*x; -b + A*x + ones(m,1)];
    [y, ~, exitflag, ~, lambda] = quadprog(Q, -(c + D*x), lhs, rhs, [], [], [], [], [], options);
    % [y, ~, exitflag, ~, lambda] = linprog(-(c + D*x), W, b - A*x, [], [], [], [], options);
    if exitflag == -3
        error('Optimal value function does not have a solution for current x');
    else
        duals = lambda.ineqlin;
    end
    
    %% (Drop due to inefficiency) Call Gurobi to solve QP %%
    % model2.A = sparse(W);
    % model2.rhs = b - A*x;
    % model2.sense = '<';
    % model2.lb = -inf(m, 1);
    % linear_term = (c + D*x) + epsilon * (D * (A' * duals));
    % Q_total = Q/2 + epsilon/2 * (D * D');
    % model2.obj = linear_term;
    % model2.Q = sparse(-Q_total);
    % model2.modelsense = 'max';
    % params.OutputFlag = 0;
    % result2 = gurobi(model2, params);
    % if ~strcmp(result2.status, 'OPTIMAL')
    %     error('QP to obtain G_oracle not solved to optimal.');
    % end
    % soln = result2.x;

    % %% Call Matlab Opt Toolbox to solve QP %%
    options.Display = 'off';
    options.Algorithm = 'active-set';  % 'interior-point-convex' (default)/'trust-region-reflective'/'active-set'
    % intialize 2nd QP with the solution of the 1st QP
    [soln, ~, exitflag] = quadprog(Q + epsilon * (D * D'), ...
                                -(c + D*x), ...
                                lhs, ...
                                rhs, ...
                                [], [], [], [], ...
                                y, ...
                                options);
    if exitflag == -3
        error('QP to obtain G_oracle not solved to optimal.');
    end

    approx_grad = 4 * norm(x)^2 * x + D' * soln - [A; -A]' * duals;
end

function grad = optval_LICQ_subgrad(x, A, b, c, D, W, Q)
    m = length(c);
    % %% (Drop due to inefficiency) Call Gurobi to solve QP %%
    % model.A = sparse(W);
    % model.rhs = b - A*x;
    % model.sense = '<';
    % model.lb = -inf(m, 1);
    % model.obj = (c + D*x);
    % model.Q = sparse(-Q/2);
    % model.modelsense = 'max';
    % params.OutputFlag = 0;
    % result = gurobi(model, params);
    % if ~strcmp(result.status, 'OPTIMAL')
    %     error('Gurobi did not return an optimal solution.');
    % end
    % duals = result.pi;
    % grad = x + D' * result.x - A' * duals;

    %% Call Matlab Opt Toolbox to solve QP %%
    options.Display = 'off';
    % options.Algorithm = 'active-set';  % 'interior-point-convex' (default)/'trust-region-reflective'/'active-set'
    
    % y0 = zeros(m,1);
    lhs = [W; -W];
    rhs = [b - A*x; -b + A*x + ones(m,1)];
    [soln, ~, exitflag, ~, lambda] = quadprog(Q, -(c + D*x), lhs, rhs, [], [], [], [], [], options);
    % [soln, ~, exitflag, ~, lambda] = linprog(-(c + D*x), W, b - A*x, [], [], [], [], options);
    if exitflag == -3
        error('Optimal value function does not have a solution for current x');
    else
        duals = lambda.ineqlin;
    end
    grad = 4 * norm(x)^2 * x + D' * soln - [A; -A]' * duals;
end