classdef Polyak
    properties
        params              % Optimization parameters (vector)
        num_oracle_iter     % Counter for oracle calls in each iteration
        loss                % Current function value
        opt_f               % Optimal function value (must be provided)
        % c1                  % Line search parameter (not used in Polyak, included for future extension)
        % c2                  % Line search parameter (not used in Polyak, included for future extension)
        stepsize            % Effective stepsize
        name
    end

    methods
        function obj = Polyak(initial_params, opt_f) %, ls_params)
            % Constructor to initialize the optimizer
            if nargin < 2 || isinf(opt_f)
                error("Polyak requires opt_value. Input opt_value!");
            end
            % if nargin < 3
            %     ls_params = struct();
            % end
            
            obj.params = initial_params;
            obj.num_oracle_iter = 0;
            obj.loss = inf;
            obj.opt_f = opt_f;

            % Set line search parameters if provided
            % obj.c1 = getOrDefault(ls_params, 'c1', 0);
            % obj.c2 = getOrDefault(ls_params, 'c2', 0.5);
            obj.stepsize = 0;
            obj.name = 'Polyak';
        end

        function obj = step(obj, func_value, subgrad)
            % Perform a single optimization step using the Polyak method

            % Evaluate function and gradient
            obj.loss = func_value(obj.params);
            d = subgrad(obj.params);

            % Update oracle call count
            obj.num_oracle_iter = obj.num_oracle_iter + 1;

            % Compute squared norm of the gradient
            nrmd2 = d' * d;

            % Compute step size if gradient norm is large enough
            if nrmd2 > 1e-20
                t = (obj.loss - obj.opt_f) / nrmd2;
            else
                t = 0;
            end

            % Update parameters
            obj.params = obj.params - t * d;
            obj.stepsize = t;
        end
    end
end

% function val = getOrDefault(structure, field, defaultVal)
%     % Helper function to get a field from a structure or return a default value
%     if isfield(structure, field)
%         val = structure.(field);
%     else
%         val = defaultVal;
%     end
% end
