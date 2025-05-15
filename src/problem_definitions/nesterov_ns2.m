%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% nesterov_2nd_bad_function.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [func_value, subgrad, G_oracle, flags, ...
            f1, f2, g1, g2] = nesterov_ns2(dimension)
% NESTEROV_2ND_BAD_FUNCTION
%   f(x) = (x(1)-1)^2/4 + sum_{i=2}^n |x(i) - 2*x(i-1)^2 + 1|
%
% Outputs:
%   func_value : function handle for f(x)
%   subgrad    : function handle nesterov_2nd_bad_subgrad(x) that returns the exact subgradient
%   G_oracle   : function handle nesterov_2nd_bad_G_oracle(x, epsilon) that approximates the subgradient
%   flags      : structure with additional information
%   --------------------------------------------------
%   optional DC oracle for f = f1 (convex) - f2 (convex)
%   f1, f2     : function handle for f1(x) & f2(x)
%   g1, g2     : subgradient oracle for f1(x) & f2(x)

    if nargin < 1, dimension = 10; end
    
    func_value = @(x) nesterov_2nd_bad_func_value(x);
    G_oracle   = @(x, epsilon) nesterov_2nd_bad_G_oracle(x, epsilon);
    subgrad    = @(x) nesterov_2nd_bad_subgrad(x);
    
    f1 = @(x) nesterov_2nd_bad_func_value(x) + 2 * norm(x(1:end-1))^2;
    f2 = @(x) 2 * norm(x(1:end-1))^2;
    g1 = @(x) nesterov_2nd_bad_subgrad(x) + 4 * [x(1:end-1); 0];
    g2 = @(x) 4 * [x(1:end-1); 0];

    flags.function_name = 'nesterov_ns2';
    flags.opt_val = 0;
    flags.opt_soln = ones(dimension, 1);
    flags.nb_parameters = dimension;
    flags.params_for_legend = [];
end

function val = nesterov_2nd_bad_func_value(x)
    term1 = (x(1) - 1)^2 / 4;
    if length(x) > 1
        term2 = sum(abs(x(2:end) - 2 * x(1:end-1).^2 + 1));
    else
        term2 = 0;
    end
    val = term1 + term2;
end

function approx_grad = nesterov_2nd_bad_G_oracle(x, epsilon)
    n = length(x);
    const_grad = zeros(n,1);
    if n == 1
        const_grad(1) = 0.5 * (x(1) - 1);
    elseif n == 2
        const_grad(1) = 4.5 * x(1) - 0.5;
        const_grad(2) = -1;
    else
        const_grad(1) = 4.5 * x(1) - 0.5;
        const_grad(2:n-1) = -1 + 4 * x(2:n-1);
        const_grad(n) = -1;
    end
    A = zeros(n-1, n);
    for i = 1:n-1
        A(i,i) = -8 * x(i);
        A(i,i+1) = 2;
    end
    c_vec = zeros(n-1, 1);
    for i = 1:n-1
        c_vec(i) = 2 * (x(i+1) - 2 * x(i)^2 + 1) - epsilon * (A(i,:) * const_grad);
    end
    % %% (Drop due to inefficiency) Call Gurobi to solve QP with box constraints %%
    % model.A = sparse(eye(n-1));
    % model.rhs = ones(n-1,1);
    % model.sense = '<';
    % model.lb = zeros(n-1,1);
    % % model.ub = ones(n-1,1);
    % model.obj = -c_vec;
    % model.Q = sparse(epsilon/2 * (A * A')); 
    % params.OutputFlag = 0;
    % result = gurobi(model, params);
    % soln = result.x;
    
    %% Call Matlab Opt Toolbox to solve QP with box constraints %%
    options.Display = 'off';
    options.Algorithm = 'active-set';  % 'interior-point-convex' (default)/'trust-region-reflective'/'active-set'
    % For a 'trust-region-reflective' equality-constrained problem, the default value is 1e-6.
    % For a 'trust-region-reflective' bound-constrained problem, the default value is 100*eps, about 2.2204e-14.
    % For the 'interior-point-convex' and 'active-set' algorithms, the default value is 1e-8.
    % options.OptimalityTolerance = 1e-14; 
    
    % Set initial point based on active index
    y0 = zeros(n-1,1);
    for i=1:n-1
        y0(i) = (sign(1+x(i+1)-2*x(i)^2) + 1)/2;
    end
    soln = quadprog(epsilon * (A * A'), -c_vec, [], [], [], [], zeros(n-1,1), ones(n-1,1), y0, options);

    %% Update regularized-subgradient %%
    approx_grad = const_grad + (soln' * A)';
end

function grad = nesterov_2nd_bad_subgrad(x)
    % n = length(x);
    % const_grad = zeros(n,1);
    % if n == 1
    %     const_grad(1) = 0.5 * (x(1) - 1);
    % elseif n == 2
    %     const_grad(1) = 4.5 * x(1) - 0.5;
    %     const_grad(2) = -1;
    % else
    %     const_grad(1) = 4.5 * x(1) - 0.5;
    %     const_grad(2:n-1) = -1 + 4 * x(2:n-1);
    %     const_grad(n) = -1;
    % end
    % A = zeros(n-1, n);
    % for i = 1:n-1
    %     A(i,i) = -8 * x(i);
    %     A(i,i+1) = 2;
    % end
    % sol = solve_nesterov_qp(A, n);
    % grad = const_grad + (sol' * A)';
    n = length(x);
    grad = zeros(n, 1);
    grad(1) = -(1-x(1))/2;
    for i=1:n-1
        r = sign(1+x(i+1)-2*x(i)^2);
        grad(i+1) = grad(i+1) + r;
        grad(i) = grad(i) - 4*x(i)*r;
    end
end

