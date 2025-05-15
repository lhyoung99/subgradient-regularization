%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% min_of_smooth_random.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [func_value, subgrad, G_oracle, flags, ...
            f1, f2, g1, g2] = min_of_smooth_random(nb_functions, nb_features, dimension, seed)
% MIN_OF_SMOOTH_RANDOM
%   f(x) = min_{i in [m]} { || A_i x - b_i ||^2 }
%        = min_{i in [m]} {  ||b_i||^2 +  g_i * x + 0.5 * x' * H_i * x }
% Outputs:
%   func_value : function handle for f(x)
%   subgrad    : function handle min_of_smooth_random_subgrad(x) that returns the exact subgradient
%   G_oracle   : function handle min_of_smooth_random_G_oracle(x, epsilon) that approximates the subgradient
%   flags      : structure with additional information

    rng(seed, "twister");

    % Generate matrix A_i with linearly independent columns
    A = zeros(nb_features, dimension, nb_functions);
    for i = 1:nb_functions
        tmp = randn(nb_features, dimension) / sqrt(nb_features);
        [Q, R] = qr(tmp, 0);
        % Check each diagonal entry of R and perturb if too small
        for j = 1:dimension
            if abs(R(j,j)) < eps
                R(j,j) = 1e-6;  % Small perturbation value
                % If sign returns 0, you may simply set R(j,j) = epsilon
            end
        end
        A(:,:,i) = Q * R;
    end
    
    soln = randn(dimension, 1) / sqrt(dimension);
    b = squeeze(pagemtimes(A, soln)); % (nb_features x nb_functions)

    H = zeros(dimension, dimension, nb_functions);
    g = zeros(dimension, nb_functions);
    for i = 1:nb_functions
        T = squeeze(A(:,:,i));
        H(:,:,i) = T' * T;
        g(:,i) = -T' * b(:,i);
    end

    func_value = @(x) min_of_smooth_random_func_value(x, A, b);
    G_oracle   = @(x, epsilon) min_of_smooth_random_G_oracle(x, epsilon, H, g, b);
    subgrad    = @(x) min_of_smooth_random_subgrad(x, A, b);
    
    f1 = func_value;
    f2 = @(x) 0;
    g1 = subgrad;
    g2 = @(x) zeros(dimension, 1);
    
    flags.name = '$f(x) = \min_{i\in[m]}\{\frac{1}{2}||A_i x - b_i||\}$';
    flags.function_name = 'min_of_smooth_random';
    flags.nb_parameters = dimension;
    flags.params_for_legend = [];
    flags.opt_val = 0;
    flags.opt_soln = soln;
end

function val = min_of_smooth_random_func_value(x, A, b)
    residual = squeeze(pagemtimes(A, x)) - b;
    val = 0.5 * min(sum(residual.^2, 1));
end

function approx_grad = min_of_smooth_random_G_oracle(x, epsilon, H, g, b)
    nb_functions = size(H, 3);
    
    Hx = squeeze(pagemtimes(H, x)); 
    Q = g' + Hx';
    c_vec = 0.5 * sum(b.^2, 1)' + g' * x + 0.5 * Hx' * x;
    
    %% Call Matlab Opt Toolbox to solve QP %%
    options.Display = 'off';
    options.Algorithm = 'active-set';  % 'interior-point-convex' (default)/'trust-region-reflective'/'active-set'
    % For a 'trust-region-reflective' equality-constrained problem, the default value is 1e-6.
    % For a 'trust-region-reflective' bound-constrained problem, the default value is 100*eps, about 2.2204e-14.
    % For the 'interior-point-convex' and 'active-set' algorithms, the default value is 1e-8.
    % options.OptimalityTolerance = 1e-14; 

    y0 = ones(nb_functions,1) / nb_functions;
    y_soln = quadprog(epsilon * (Q * Q'), c_vec, [], [], ones(1, nb_functions), 1, zeros(nb_functions,1), ones(nb_functions,1), y0, options);

    approx_grad = (y_soln' * Q)';
end

function grad = min_of_smooth_random_subgrad(x, A, b)
    residual = squeeze(pagemtimes(A, x)) - b;
    vals = sum(residual.^2, 1);
    [~, idx] = min(vals);
    % fprintf("active index: %d\n",idx);
    A_idx = squeeze(A(:, :,idx));
    grad = A_idx' * (A_idx * x - b(:,idx));
end