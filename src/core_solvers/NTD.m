classdef NTD
    properties
        params              % Optimization parameters (vector)
        num_oracle          % Number of oracle calls
        loss                % Current function value
        stationary_measure  % Measure of stationarity
        sigma_increase      % Whether to increase sigma adaptively
        use_trust_region    % Whether to use trust region
        s                  % Scale factor for trust region
        s_scale_factor      % Scaling factor for s
        nb_increasing_steps % Max allowable increasing steps
        % sigma              % Step size
        stepsize            % Effective stepsize
        opt_f              % Optimal function value if known
        name
    end

    methods
        function obj = NTD(initial_params, opt_f, method_params)
            % Constructor to initialize optimizer 
            obj.params = initial_params;
            obj.num_oracle = [];
            obj.loss = inf;
            obj.stationary_measure = [];
            obj.sigma_increase = method_params.adaptive_grid_size;
            obj.use_trust_region = method_params.use_trust_region;
            obj.s_scale_factor = method_params.s_scale_factor;
            obj.nb_increasing_steps = inf;
            obj.s = [];
            obj.opt_f = opt_f;
            % obj.sigma = [];
            obj.stepsize = 0;
            obj.name = 'NTD';
        end

        function [loss, grad] = directional_evaluate(~, func_value, subgrad, x, t, d)
            % Evaluate function and subgradient in a given direction
            % obj.params = x + t * d;
            x_new = x + t * d;
            loss = func_value(x_new);
            grad = subgrad(x_new);
        end

        function g = optimal_average(~, g, hatg)
            % Compute an optimal average between two gradients
            y = hatg - g;
            dp = g' * y;
            nrmsquare = y' * y;
            if nrmsquare == 0
                g = g;
            else
                weight = max(min(-dp / nrmsquare, 1), 0);
                g = (1 - weight) * g + weight * hatg;
            end
        end

        function [g, nb1_oracle, min_sub_norm] = NDescent(obj, func_value, subgrad, x, g, loss, sigma, T)
            % Nonsmooth descent method
            best_g = g;
            min_value = loss;
            nb1_oracle = 0;
            for i = 1:T
                nrmg = norm(g);
                if nrmg <= 1e-20, break; end
                [f_new, ~] = obj.directional_evaluate(func_value, subgrad, x, -sigma / nrmg, g);
                nb1_oracle = nb1_oracle + 1;
                if f_new < min_value
                    best_g = g;
                    min_value = f_new;
                end
                if f_new <= loss - sigma * nrmg / 8
                    best_g = g;
                    break;
                end
                t = rand();
                [~, hatg] = obj.directional_evaluate(func_value, subgrad, x, -t * sigma / nrmg, g);
                nb1_oracle = nb1_oracle + 1;
                g = obj.optimal_average(g, hatg);
            end
            min_sub_norm = norm(g);
        end

        function [g, nb2_oracle, min_sub_norm] = TDescent(obj, func_value, subgrad, x, g, loss, sigma, T)
            % Trust-region descent method
            best_g = g;
            min_value = loss;
            nb2_oracle = 0;
            for i = 1:T
                nrmg = norm(g);
                if nrmg <= 1e-20, break; end
                [f_new, hatg] = obj.directional_evaluate(func_value, subgrad, x, -sigma / nrmg, g);
                nb2_oracle = nb2_oracle + 1;
                if f_new < min_value
                    best_g = g;
                    min_value = f_new;
                end
                if f_new <= loss - sigma * nrmg / 8
                    best_g = g;
                    break;
                end
                g = obj.optimal_average(g, hatg);
            end
            min_sub_norm = norm(g);
        end

        function obj = step(obj, func_value, subgrad)
            % Perform an optimization step
            
            x = obj.params;
            obj.loss = func_value(x);
            grad = subgrad(x);

            % Line search initialization
            K = length(obj.num_oracle); % Iteration counter
            T = K;
            G = min(K, 50); % Set upper limit to G

            % Run line search
            if K == 0
                obj.s = norm(grad);
            end
            [a, sigma, nb_oracle, R_k] = obj.linesearch(func_value, subgrad, x, grad, obj.loss, G, T);
            % obj.sigma = sigma;

            % Update parameters
            nrm_a = norm(a);
            if nrm_a >= eps
                t = sigma / nrm_a;
                obj.params = obj.params - t * a;
                obj.stepsize = t;
            end

            % Store iteration statistics
            obj.num_oracle = [obj.num_oracle, nb_oracle];
            obj.stationary_measure = [obj.stationary_measure, R_k];
        end

        function [best_g, best_sigma, nb_oracle, R_k] = linesearch(obj, func_value, subgrad, x, g, loss, G, T)
            % Line search strategy for step size selection
            
            min_value = loss;
            best_g = g;
            best_hatg = g;
            best_sigma = 0;
            v = g; 
            nb_oracle = 0;
            constant = max(norm(g), obj.s * obj.s_scale_factor);
            
            best_idx_so_far = 0;
            if isempty(obj.opt_f) || obj.opt_f == inf
                dist_est = 1;
            else
                dist_est = 10 * ((loss - obj.opt_f) / constant);
            end
            nrmv = norm(v);

            % if R_k is overflow, then we set R_k = 1e10
            if nrmv > 60
                R_k = inf;
            else
                R_k = 2^(nrmv);
            end

            sigma_multiplier = 1;

            sigma = dist_est * sigma_multiplier * 2^(-G);
            nb_increasing_steps = 0;
            min_sub_norm = norm(g);
            i = 0;

            while i < G && (sigma_multiplier * min_sub_norm / constant >= sigma || ~obj.use_trust_region)
                i = i + 1;
                sigma = dist_est * sigma_multiplier * 2^(-(G - i));
                [u, no1, min_sub_norm] = obj.NDescent(func_value, subgrad, x, v, loss, sigma, T);
                nb_oracle = nb_oracle + no1;
                [v, no2, min_sub_norm] = obj.TDescent(func_value, subgrad, x, u, loss, sigma, T);
                nb_oracle = nb_oracle + no2;
                nrmv = norm(v);
                % Take the maximum of np.power(nrmv, 2.0) and sigma
                max_of_vals = max(min_sub_norm^2, sigma);
                % Update R_k to be the minimum of R_k and max_of_vals
                R_k = min(R_k, max_of_vals);
                if nrmv < 1e-20, break; end

                [f_new, hatg] = obj.directional_evaluate(func_value, subgrad, x, -sigma / nrmv, v);
                nb_oracle = nb_oracle + 1;
                if nrmv < 1e-20, break; end
                if f_new < min_value
                    min_value = f_new;
                    best_g = v;
                    best_hatg = hatg;
                    best_sigma = sigma;
                    best_idx_so_far = i;
                else
                    if min_value < loss
                        nb_increasing_steps = nb_increasing_steps + 1;
                        if nb_increasing_steps > obj.nb_increasing_steps
                            break;
                        end
                    end
                end

                if obj.sigma_increase
                    if best_idx_so_far == G
                        G = G+1;
                        sigma_multiplier = sigma_multiplier * 10;
                    end
                end
            end

            nrm_best_g = norm(best_g); 
        end
    end
end
