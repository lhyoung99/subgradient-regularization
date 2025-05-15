function results = SRD_adapt(...
    x0, ... 
    v, ...
    G_func, ...
    flags, ...
    method_params ...
)
%{
Subgradient-Regularized Descent (Adaptive Version) with Oracle Call Limitation and Information Display.

Parameters:
    x0: Initial point (vector). 
    v: Objective function handle (minimization). For example, @(x) ...
    G_func: Function handle that takes (x, epsilon) and returns the descent direction.
    flags: Struct containing settings and flags for optimization. 
    method_params: Struct containing parameter setting for algorithms

Returns:
    A struct containing optimization results and statistics.
%}

epsilon_0_0 = method_params.epsilon_0_0;
epsilon_tol = method_params.epsilon_tol;
nu_tol = method_params.nu_tol;
nu_0 = method_params.nu_0;
alpha_nu = method_params.alpha_nu; % Reduction factor (0,1)
alpha_epsilon = method_params.alpha_epsilon; % Reduction factor (0,1)
beta = method_params.beta; % Armijo line search param (0,1)
obj_tol = method_params.obj_tol;
max_oracle_call = method_params.max_oracle_call;
max_time_seconds = method_params.max_time_seconds;
print_frequency = method_params.print_frequency;

dimension = length(x0);
x_k = x0;
epsilon_k_0 = epsilon_0_0;
epsilon_test = 1;   % Adaptive subgradient parameter
eta_k = NaN;        % Will be found in line search
nu_k = nu_0;
num_oracle = 0;
num_oracle_list = [];
objective_values_list = [];
distance_list = [];
history = [];
stepsize = 0;
stepsize_list = [];
header_printed = false;

method_name = 'SRDescent-adapt';

opt_val = [];
if isfield(flags, 'opt_val')
    opt_val = flags.opt_val;
end

opt_soln = [];
if isfield(flags, 'opt_soln')
    opt_soln = flags.opt_soln;
end

start_time = tic;
elapsed_time_list = [];

k = 0;
t = 0;
obj_current = v(x_k);

if ~isempty(opt_soln)
    distance_current = norm(x_k - opt_soln);
else
    distance_current = [];
end

num_oracle_list(end+1) = num_oracle;
objective_values_list(end+1) = obj_current; 
history = [history, x_k(:)];
stepsize_list = [stepsize_list, stepsize];
distance_list(end+1) = distance_or_empty(distance_current);
elapsed_time_list(end+1) = toc(start_time);

while true
    i = 0;
    while true
        % epsilon_{k,i} = epsilon_{k,0} * 2^{-i}
        epsilon_k_i = epsilon_k_0 * (2^(-i));
        
        % Compute descent direction g_k_i
        g_k_i = G_func(x_k, epsilon_k_i);
        num_oracle = num_oracle + 1;
        
        % Record data
        % num_oracle_list(end+1) = num_oracle;
        % objective_values_list(end+1) = obj_current;
        % history = [history, x_k(:)];
        % distance_list(end+1) = distance_or_empty(distance_current);
        % elapsed_time_list(end+1) = toc(start_time);
        
        % Check for termination
        if should_terminate(...
                elapsed_time_list(end), ...
                max_time_seconds, ...
                num_oracle, ...
                max_oracle_call, ...
                obj_current, ...
                opt_val, ...
                obj_tol, ...
                epsilon_tol, ...
                nu_tol, ...
                epsilon_k_i, ...
                g_k_i)
            
            update_and_log_progress(...
                k, i, header_printed, g_k_i, [], epsilon_k_0, [], [], ...
                num_oracle, obj_current, distance_current);
            
            results = compile_statistics(method_name, x_k, num_oracle_list, ...
                elapsed_time_list, objective_values_list, distance_list, history, flags);
            return;
        end
        
        % Perform line search
        eta_0_local = epsilon_k_0; 
        [eta_candidate, nb_oracles, x_k_new, obj_current, ls_stop_flag] = perform_line_search(...
            i, eta_0_local, beta, x_k, g_k_i, v, obj_current);
        
        num_oracle = num_oracle + nb_oracles;
        
        % Check if line search "failed" below machine precision
        if ls_stop_flag
            update_and_log_progress(...
                k, i, header_printed, g_k_i, [], epsilon_k_0, [], [], ...
                num_oracle, obj_current, distance_current);
            
            % Record data
            num_oracle_list(end+1) = num_oracle;
            objective_values_list(end+1) = obj_current;
            history = [history, x_k(:)];
            stepsize_list = [stepsize_list, stepsize];
            distance_list(end+1) = distance_or_empty(distance_current);
            elapsed_time_list(end+1) = toc(start_time);
            results = compile_statistics(method_name, x_k, num_oracle_list, ...
                elapsed_time_list, objective_values_list, distance_list, history, flags);
            return;
        end

        if ~isempty(x_k_new)
            % Successful step
            stepsize = eta_candidate * norm(g_k_i);
            if ~isempty(opt_soln)
                distance_current = norm(x_k_new - opt_soln);
            end
            
            % If norm(g_k_i) <= nu_k, then do additional checks with epsilon_test
            if norm(g_k_i) <= nu_k
                t = t + 1;
                epsilon_test = t^(-0.25);
                g_test = G_func(x_k, epsilon_test);
                num_oracle = num_oracle + 1;
                
                [k, header_printed] = update_and_log_progress(...
                    k, i, header_printed, g_k_i, g_test, epsilon_k_0, ...
                    epsilon_test, eta_candidate, num_oracle, ...
                    obj_current, distance_current);
                
                % Additional stationarity check
                if ~isempty(epsilon_tol) && ~isempty(nu_tol) && (epsilon_test <= epsilon_tol) && (norm(g_test) <= nu_tol)
                    disp('Approximate Clarke stationarity reached (Testing subgrad).');
                    results = compile_statistics(method_name, x_k, num_oracle_list, ...
                        elapsed_time_list, objective_values_list, distance_list, history, flags);
                    return;
                end
                
                nu_k = alpha_nu * nu_k;

                % Adjust epsilon_k_0 based on the Ratio Test
                lhs = epsilon_test * norm(g_test);
                rhs = (1 / epsilon_k_0) * sqrt(epsilon_k_i * norm(g_k_i));
                if lhs >= rhs
                    epsilon_k_0 = alpha_epsilon * epsilon_k_0;
                end
            else
                % If norm(g_k_i) > nu_k, no g_test needed
                g_test = [];
                if mod(k, print_frequency) == 0
                    [k, header_printed] = update_and_log_progress(...
                        k, i, header_printed, g_k_i, g_test, ...
                        epsilon_k_0, epsilon_test, eta_candidate, ...
                        num_oracle, obj_current, distance_current);
                else
                    k = k + 1;
                end
            end
            
            x_k = x_k_new;
            
            % Record data
            num_oracle_list(end+1) = num_oracle;
            objective_values_list(end+1) = obj_current;
            history = [history, x_k_new(:)];
            stepsize_list = [stepsize_list, stepsize];
            elapsed_time_list(end+1) = toc(start_time);
            distance_list(end+1) = distance_or_empty(distance_current);
            
            break; % Exit the inner loop
        end
        
        i = i + 1; % Try smaller epsilon_k_i
    end
end
end % end of subgradient_regularized_descent_adapt


%% ===================== Helper functions =====================
function val = distance_or_empty(distance_current)
    if isempty(distance_current)
        val = NaN;
    else
        val = distance_current;
    end
end

function stop_flag = should_terminate(...
    elapsed_time, ...
    max_time_seconds, ...
    num_oracle, ...
    max_oracle_call, ...
    obj_current, ...
    opt_val, ...
    obj_tol, ...
    epsilon_tol, ...
    nu_tol, ...
    epsilon_k_i, ...
    g_k_i ...
)
    stop_flag = false;
    if (elapsed_time > max_time_seconds) || (num_oracle >= max_oracle_call)
        disp('Maximum oracle calls or time limit reached.');
        stop_flag = true;
        return;
    end
    
    if ~isempty(opt_val) && ~isempty(obj_tol) && ((obj_current - opt_val) <= obj_tol)
        disp('Objective tolerance reached.');
        stop_flag = true;
        return;
    end
    
    if ~isempty(epsilon_tol) && ~isempty(nu_tol) && ~isempty(g_k_i) && (epsilon_k_i <= epsilon_tol) && (norm(g_k_i) <= nu_tol)
        disp('Approximate Clarke stationarity reached.');
        stop_flag = true;
        return;
    end

    % if epsilon_k_i <= eps
    %     disp('Terminate because epsilon is decreased under the machine precision.');
    %     stop_flag = true;
    %     return;
    % end
end

function [eta, nb_oracles, x_k_new, obj_current, ls_stop_flag] = perform_line_search(i, eta_0_local, beta, x_k, g_k_i, v, obj_current)
    eta = [];
    x_k_new = [];
    obj_best = inf;
    nb_oracles = 0;
    ls_stop_flag = false;
    eta_lowbound = eps * norm(x_k)/norm(g_k_i);
    for j = 0:i
        eta_candidate = eta_0_local * (2^(-j));
        x_trial = x_k - eta_candidate*g_k_i;
        obj_new = v(x_trial);
        nb_oracles = nb_oracles + 1;

        if obj_new < obj_best && obj_new <= obj_current - eta_candidate * beta * (norm(g_k_i)^2)
            % % Stop if improvement is below machine precision (eps ~ 2e-16)
            % tol = eps * max(abs(obj_current), abs(obj_new));
            % if (0 <= obj_current - obj_new) && (obj_current - obj_new < tol)
            %     disp('Terminate: The improvement of line search is below machine precision.');
            %     ls_stop_flag = true;
            %     return;
            % end

            eta = eta_candidate;
            x_k_new = x_trial;
            obj_current = obj_new;
            obj_best = obj_new;
        elseif eta_candidate < eta_lowbound
            % Stop if improvement is below machine precision (eps ~ 2e-16)
            disp('Terminate: The improvement of line search is below machine precision.');
            ls_stop_flag = true;
            return;
        end
    end
end

function [k_new, header_printed_out] = update_and_log_progress(...
    k, i, header_printed, g_k_i, g_test, epsilon_k_0, ...
    epsilon_test, eta_k, num_oracle, obj_current, distance_current ...
)
    if ~header_printed
        fprintf('%-15s%-15s%-5s%-15s%-15s%-15s%-15s%-15s%-15s%-20s%-15s\n', ...
            'k', '||x_k - x*||', 'i', 'eta', 'epsilon_k_0', 'epsilon_k_i', ...
            '||g_k_i||', 'epsilon_test', '||g_test||', 'v(x_k)', 'num_oracle');
        fprintf('%s\n', repmat('-', 1, 140));
        header_printed = true;
    end
    
    if isempty(distance_current)
        distance_str = '--';
    else
        distance_str = num2str(distance_current, '%.4g');
    end
    
    if isempty(eta_k)
        eta_str = '--';
    else
        eta_str = sprintf('%.4e', eta_k);
    end
    
    if isempty(epsilon_test)
        eps_test_str = '--';
    else
        eps_test_str = sprintf('%.4e', epsilon_test);
    end
    
    if isempty(g_test)
        g_test_norm_str = '--';
    else
        g_test_norm_str = sprintf('%.4e', norm(g_test));
    end
    
    epsilon_k_i_val = epsilon_k_0*(2^(-i));
    
    fprintf('%-15d%-15s%-5d%-15s%-15.4e%-15.4e%-15.4e%-15s%-15s%-20.4e%-15d\n', ...
        k+1, distance_str, i+1, eta_str, epsilon_k_0, ...
        epsilon_k_i_val, norm(g_k_i), eps_test_str, ...
        g_test_norm_str, obj_current, num_oracle);
    
    k_new = k + 1;
    header_printed_out = header_printed;
end

function results_struct = compile_statistics(...
    method_name, x_k, num_oracle_list, elapsed_time_list, ...
    objective_values_list, distance_list, history, flags ...
)
    results_struct.name_of_method = method_name;
    results_struct.x = x_k;
    results_struct.num_oracle_list = num_oracle_list;
    results_struct.elapsed_time_list = elapsed_time_list;
    results_struct.objective_values_list = objective_values_list;
    results_struct.distance_list = distance_list;
    results_struct.history = history;
    results_struct.flags = flags;
end
