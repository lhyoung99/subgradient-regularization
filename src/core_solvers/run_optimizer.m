function [optimizer, statistics] = run_optimizer(optimizer, func_value, subgrad, flags, method_params)
    
    max_oracle_call = method_params.max_oracle_call;
    max_time_seconds = method_params.max_time_seconds;
    print_frequency = method_params.print_frequency;
    save_frequency = method_params.save_frequency;
    opt_val = flags.opt_val;
    obj_tol = method_params.obj_tol;

    opt_solution = [];
    if isfield(flags, 'opt_soln')
        opt_solution = flags.opt_soln;
    end
    if isempty(opt_val) && isfield(flags, 'opt_val')
        opt_val = flags.opt_val;
    end
    % if isempty(opt_val) && ~isempty(obj_tol)
    %     error('Cannot exit based on tolerance with no opt_val. Input opt_val or set obj_tol to empty.');
    % end

    param = optimizer.params; 
    name_of_method = optimizer.name;
    numel_param = flags.nb_parameters;
    hasNumOracle = isprop(optimizer, 'num_oracle');
    hasStationaryMeasure = isprop(optimizer, 'stationary_measure');

    % History storage
    history = zeros(numel_param, 0);
    stepsize_list = [];
    loss_value = [];
    loss_best = Inf;
    distance_list = [];
    num_oracle_list = [];
    elapsed_time = [];
    
    if hasStationaryMeasure
        fprintf('+--------------+------------------+------------------+\n');
        fprintf('| Oracle Calls |   Best Obj Val   |   Stationarity   |\n');
        fprintf('+--------------+------------------+------------------+\n');
    else
        fprintf('+--------------+------------------+\n');
        fprintf('| Oracle Calls |   Best Obj Val   |\n');
        fprintf('+--------------+------------------+\n');
    end

    % Timer start
    start_time = tic;
    nb_oracles = 0;
    % stepsize = 0;
    print_counter = 0;

    while nb_oracles < max_oracle_call
        % Optimizer step update
        optimizer = optimizer.step(func_value, subgrad);
        param = optimizer.params;
        loss = optimizer.loss;
        stepsize = optimizer.stepsize;

        % Update oracle call count
        if hasNumOracle
            nb_oracles = sum(optimizer.num_oracle);
        else
            nb_oracles = nb_oracles + 1;
        end
        
        % Compute function value
        % loss = func_value(param);
        loss_best = min([loss_best, loss]); 

        % Compute distance to optimal solution if applicable
        if ~isempty(opt_solution)
            distance = norm(param - opt_solution);
        end
        
        % Record elapsed time
        current_time = toc(start_time);

        % Store history
        if mod(print_counter, save_frequency) == 0
            loss_value = [loss_value, loss_best];
            history = [history, param];
            stepsize_list = [stepsize_list, stepsize];
            elapsed_time = [elapsed_time, current_time];
            num_oracle_list = [num_oracle_list, nb_oracles];
            if ~isempty(opt_solution)
                distance_list = [distance_list, distance];
            end
        end

        % Print progress if required
        if mod(print_counter, print_frequency) == 0
            if hasStationaryMeasure
                fprintf('| %-12d | %-16.4e | %-16.4e |\n', ...
                        nb_oracles, loss_best, optimizer.stationary_measure(end));
            else
                fprintf('| %-12d | %-16.4e |\n', ...
                        nb_oracles, loss_best);
            end
        end

        % Stopping conditions
        if ~isempty(opt_val) && ~isempty(obj_tol) && (loss_best - opt_val < obj_tol)
            fprintf('Objective tolerance reached.\n');
            x_result = param;
            % Store history & Print
            if mod(print_counter, save_frequency) ~= 0
                loss_value = [loss_value, loss_best];
                history = [history, param];
                stepsize_list = [stepsize_list, stepsize];
                elapsed_time = [elapsed_time, current_time];
                num_oracle_list = [num_oracle_list, nb_oracles];
                if ~isempty(opt_solution)
                    distance_list = [distance_list, distance];
                end
            end
            if mod(print_counter, print_frequency) ~= 0
                if hasStationaryMeasure
                    fprintf('| %-12d | %-16.4e | %-16.4e |\n', ...
                            nb_oracles, loss_best, optimizer.stationary_measure(end));
                else
                    fprintf('| %-12d | %-16.4e |\n', ...
                            nb_oracles, loss_best);
                end
            end
            break;
        elseif current_time > max_time_seconds
            fprintf('Maximum time limit reached.\n');
            x_result = param;
            % Store history & Print
            if mod(print_counter, save_frequency) ~= 0
                loss_value = [loss_value, loss_best];
                history = [history, param];
                stepsize_list = [stepsize_list, stepsize];
                elapsed_time = [elapsed_time, current_time];
                num_oracle_list = [num_oracle_list, nb_oracles];
                if ~isempty(opt_solution)
                    distance_list = [distance_list, distance];
                end
            end
            if mod(print_counter, print_frequency) ~= 0
                if hasStationaryMeasure
                    fprintf('| %-12d | %-16.4e | %-16.4e |\n', ...
                            nb_oracles, loss_best, optimizer.stationary_measure(end));
                else
                    fprintf('| %-12d | %-16.4e |\n', ...
                            nb_oracles, loss_best);
                end
            end
            break;
        elseif nb_oracles >= max_oracle_call
            fprintf('Maximum oracle calls reached.\n');
            x_result = param;
            % Store history & Print
            if mod(print_counter, save_frequency) ~= 0
                loss_value = [loss_value, loss_best];
                history = [history, param];
                stepsize_list = [stepsize_list, stepsize];
                elapsed_time = [elapsed_time, current_time];
                num_oracle_list = [num_oracle_list, nb_oracles];
                if ~isempty(opt_solution)
                    distance_list = [distance_list, distance];
                end
            end
            if mod(print_counter, print_frequency) ~= 0
                if hasStationaryMeasure
                    fprintf('| %-12d | %-16.4e | %-16.4e |\n', ...
                            nb_oracles, loss_best, optimizer.stationary_measure(end));
                else
                    fprintf('| %-12d | %-16.4e |\n', ...
                            nb_oracles, loss_best);
                end
            end
            break;
        end

        print_counter = print_counter + 1;
    end

    % Store statistics
    statistics.x = x_result;
    statistics.elapsed_time_list = elapsed_time;
    statistics.objective_values_list = loss_value;
    statistics.distance_list = distance_list;
    statistics.num_oracle_list = num_oracle_list;
    statistics.history = history;
    statistics.stepsize_list = stepsize_list;
    statistics.flags = flags;

    % Rename method for output consistency
    if strcmp(name_of_method, 'NTD')
        statistics.name_of_method = 'NTDescent';
    elseif strcmp(name_of_method, 'Polyak')
        statistics.name_of_method = 'Polyak';
    else
        statistics.name_of_method = name_of_method;
    end
end
