%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% counter_example_2D.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [func_value, subgrad, G_oracle, flags] = counter_example_2D()
% COUNTER_EXAMPLE_2D
%   f(x) = max{-100, 2*x1+3*x2, -2*x1+3*x2, 5*x1+2*x2, -5*x1+2*x2}
%
% Outputs:
%   func_value : function handle for f(x)
%   subgrad    : function handle counter_example_2D_subgrad(x) that returns the exact subgradient
%   G_oracle   : function handle counter_example_2D_G_oracle(x, epsilon) that computes an approximate subgradient
%   flags      : structure with additional information

    % (No initial x is output; the user supplies x to the returned function handles.)
    func_value = @(x) max([-100, 2*x(1)+3*x(2), -2*x(1)+3*x(2), 5*x(1)+2*x(2), -5*x(1)+2*x(2)]);
    G_oracle = @(x, epsilon) counter_example_2D_G_oracle(x, epsilon);
    subgrad = @(x) counter_example_2D_subgrad(x);
    
    flags.name = '$f(x) = \max\{-100, 2x_1+3x_2, -2x_1+3x_2, 5x_1+2x_2, -5x_1+2x_2\}$';
    flags.function_name = 'counter_example_2D';
    flags.opt_val = -100;
    flags.opt_soln = [];
    flags.dimension = 2;
    flags.nb_parameters = 2;
    flags.params_for_legend = [];
end

function approx_grad = counter_example_2D_G_oracle(x, epsilon)
    grad1    = [2, 3];
    gradNeg1 = [-2, 3];
    grad2    = [5, 2];
    gradNeg2 = [-5, 2];
    A = [0, 0; grad1; gradNeg1; grad2; gradNeg2];  % 5-by-2 matrix
    c_vec = [-100; 2*x(1)+3*x(2); -2*x(1)+3*x(2); 5*x(1)+2*x(2); -5*x(1)+2*x(2)];
    
    % QP: maximize c_vec'*y - (epsilon/2)*y'*(A*A')*y subject to sum(y)==1, y>=0.
    model.A = ones(1,5);
    model.rhs = 1;
    model.sense = '=';
    model.lb = zeros(5,1);
    model.obj = c_vec;
    model.Q = -epsilon * (A * A');
    model.modelsense = 'max';
    params.OutputFlag = 0;
    result = gurobi(model, params);
    soln = result.x;
    approx_grad = (soln' * A)';  % return as column vector
end

function grad = counter_example_2D_subgrad(x)
    vals = [-100, 2*x(1)+3*x(2), -2*x(1)+3*x(2), 5*x(1)+2*x(2), -5*x(1)+2*x(2)];
    [~, idx] = max(vals);
    switch idx
        case 1, grad = [0; 0];
        case 2, grad = [2; 3];
        case 3, grad = [-2; 3];
        case 4, grad = [5; 2];
        case 5, grad = [-5; 2];
        otherwise, grad = [0; 0];
    end
end