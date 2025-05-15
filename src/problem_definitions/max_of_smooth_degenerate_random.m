%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% max_of_smooth_random.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [func_value, subgrad, G_oracle, flags, ...
            f1, f2, g1, g2] = max_of_smooth_degenerate_random(nb_functions, dimension, seed)
% MAX_OF_SMOOTH_RANDOM
%   f(x) = max_{i in [m]} { g_i * x + 0.5 * x' * H_i * x }
%
% Outputs:
%   func_value : function handle for f(x)
%   subgrad    : function handle max_of_smooth_sharp_random_subgrad(x) that returns the exact subgradient
%   G_oracle   : function handle max_of_smooth_sharp_random_G_oracle(x, epsilon) that approximates the subgradient
%   flags      : structure with additional information

    rng(seed, "twister");
    
    % Generate weights (lambda values) and ensure they sum up to 1
    lam = rand(1, nb_functions/2);
    lam = lam / sum(lam);
    alpha = rand(1, nb_functions -1);
    alpha(end) = -sum(alpha(1:end-1));

    % Generate random vectors for g
    g_temp = randn(nb_functions - 2, dimension) / sqrt(dimension);
    gk = - (lam(1:end-1) * g_temp(1:nb_functions/2 - 1,:)) / lam(end);
    gk_affine = - (alpha(1:end-1) * g_temp) / alpha(end);
    g = [g_temp; gk; gk_affine];  % (nb_functions x dimension)

    H = zeros(dimension, dimension, nb_functions);
    tmp = randn(dimension, dimension, nb_functions) / sqrt(dimension);
    for i = 1:nb_functions
        T = squeeze(tmp(:,:,i));
        H(:,:,i) = T' * T;
    end
    
    func_value = @(x) max_of_smooth_sharp_random_func_value(x, g, H);
    G_oracle   = @(x, epsilon) max_of_smooth_sharp_random_G_oracle(x, epsilon, g, H);
    subgrad    = @(x) max_of_smooth_sharp_random_subgrad(x, g, H);
    
    f1 = func_value;
    f2 = @(x) 0;
    g1 = subgrad;
    g2 = @(x) zeros(dimension, 1);

    flags.name = '$f(x) = \max_{i\in[m]}\{\frac{1}{2}x^\top A_i x + b_i^\top x\}$'; 
    flags.function_name = 'max_of_smooth_degenerate_random';
    flags.opt_val = 0;
    flags.opt_soln = zeros(dimension,1);
    flags.nb_parameters = dimension;
    flags.params_for_legend = [];
end

function val = max_of_smooth_sharp_random_func_value(x, g, H)
    Hx = squeeze(pagemtimes(H, x));
    xHx = x' * Hx;
    vals = g * x + 0.5 * xHx';
    val = max(vals);
end


function approx_grad = max_of_smooth_sharp_random_G_oracle(x, epsilon, g, H)
    % A = zeros(nb_functions, length(x));
    % c_vec = zeros(nb_functions,1);
    % for i = 1:nb_functions
    %     Hi = squeeze(H(i,:,:));
    %     A(i,:) = (Hi * x)' + g(i,:);
    %     c_vec(i) = g(i,:) * x + 0.5 * x' * (Hi * x);
    % end
    nb_functions = size(g, 1);

    Hx = squeeze(pagemtimes(H, x)); 
    A = g + Hx';
    c_vec = g * x + 0.5 * Hx' * x;

    % model.A = sparse(ones(1, nb_functions));
    % model.rhs = 1;
    % model.sense = '=';
    % model.lb = zeros(nb_functions,1);
    % model.obj = c_vec;
    % model.Q = sparse(-epsilon/2 * (A * A'));
    % model.modelsense = 'max';
    % params.OutputFlag = 0;
    % result = gurobi(model, params);
    % soln = result.x;
    % approx_grad = (soln' * A)';
    %% Call Matlab Opt Toolbox to solve QP %%
    options.Display = 'off';
    options.Algorithm = 'active-set';  % 'interior-point-convex' (default)/'trust-region-reflective'/'active-set'

    y0 = ones(nb_functions,1)/nb_functions;
    soln = quadprog(epsilon * (A * A'), -c_vec, [], [], ones(1, nb_functions), 1, zeros(nb_functions,1), ones(nb_functions,1), y0, options);

    approx_grad = (soln' * A)';
end



% function grad = max_of_smooth_random_subgrad(x, g, H, nb_functions)
%     vals = zeros(nb_functions,1);
%     grads = zeros(nb_functions, length(x));
%     for i = 1:nb_functions
%         Hi = squeeze(H(i,:,:));
%         vals(i) = g(i,:) * x + 0.5 * x' * (Hi * x);
%         grads(i,:) = g(i,:) + (Hi * x)';
%     end
%     [~, idx] = max(vals);
%     grad = grads(idx,:)';
% end

function grad = max_of_smooth_sharp_random_subgrad(x, g, H)
    % Compute all function values at once
    Hx = squeeze(pagemtimes(H, x)); 
    xHx = x' * Hx;
    vals = g * x + 0.5 * xHx';  % column vector

    % Compute all gradients at once (without transposing)
    grads = g + Hx'; 

    [~, idx] = max(vals);
    % fprintf("active index: %d\n",idx);

    % Return the corresponding gradient
    grad = grads(idx, :)';
end




