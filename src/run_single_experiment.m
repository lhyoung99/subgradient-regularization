function results_summary = run_single_experiment(example, num_initials, options_from_caller)
% RUN_SINGLE_EXPERIMENT - Executes a set of algorithms for a given problem
%
% Args:
%   example (struct): Defines the problem (e.g., example.name, example.dimension).
%   num_initials (scalar): Number of random initial points to test.
%   options_from_caller (struct): Options provided by the calling experiment script.
%                                  These will be merged with defaults from set_options.m.
%
% Returns:
%   results_summary (struct): Contains aggregated results, paths to tables if generated.


    %% 1. Finalize Options
    % Use default options if not provided
    if nargin < 3
        options_from_caller = struct();
    end
    options = set_options(options_from_caller);

    %% 2. Define Optimization Methods to Test
    available_methods = struct();
    available_methods.SRD = @SRD; 
    available_methods.SRD_adapt = @SRD_adapt; 
    available_methods.Polyak = @Polyak;
    available_methods.GradSamp = @GradSamp;
    available_methods.NTD = @NTD;
    available_methods.PBMDC = @PBMDC; 
    available_methods.BFGS = @BFGS;

    if isfield(options, 'methods_to_run') && iscell(options.methods_to_run) && ~isempty(options.methods_to_run)
        methods_to_run_names = options.methods_to_run;
        selected_methods = struct();
        for i = 1:length(methods_to_run_names)
            m_name = methods_to_run_names{i};
            if isfield(available_methods, m_name)
                selected_methods.(m_name) = available_methods.(m_name);
            else
                warning('Method "%s" requested in options.methods_to_run but not available. Skipping.', m_name);
            end
        end
        if isempty(fieldnames(selected_methods))
            error('No valid methods selected to run. Check options.methods_to_run.');
        end
    else
        selected_methods = available_methods; % Default to all available methods
    end
    method_names = fieldnames(selected_methods);
    
    %% 3. Get Problem Definition
    
    problem_func_name = example.name;
    if ~exist(problem_func_name, 'file')
        error('Problem definition function "%s.m" not found. Ensure it is in src/problem_definitions/ and on the path.', problem_func_name);
    end

    % Expected output: [func_value, subgrad, G_oracle, flags, problem_dimension, varargout for DC components]
    [problem_handles, problem_details] = local_get_problem_handles_and_details(problem_func_name, example, options);
    dimension = problem_details.dimension;
    
    %% 4. Initialize Results Storage
    % `aggregated_results` stores summaries (mean/median time, obj, success rate).
    % `all_individual_run_stats` stores the detailed 'stats' struct from every single run.
    aggregated_results_template = struct('obj_vals', [], 'times', [], 'oracle_calls', [], 'success_count', 0, 'failure_count', 0);

    aggregated_results = struct();
    aggregated_results.problem_name = example.name;
    aggregated_results.dimension = dimension;
    for k = 1:length(method_names)
        aggregated_results.(method_names{k}) = aggregated_results_template;
    end
    all_individual_run_stats = {}; % Cell array to store 'stats' from every run

    %% 5. Main Experiment Loop
    fprintf('\n%s\n', repmat('=', 1, 70));
    fprintf('===== Starting Experiment Set: %s (Dimension: %d) =====\n', example.name, dimension);
    fprintf('%s\n', repmat('=', 1, 70));
    for seed_idx = 1:num_initials
        % Generate/Set Initial Point x0
        x0_is_fixed = false;
        if isfield(options, 'x0') && length(options.x0) == dimension
            num_initials = 1;
            x0 = options.x0;
            x0_is_fixed = true;
        elseif isfield(options, 'x0')
            error('Input "options.x0" must be of length %d', dimension);
        else
            % Generate random initial points
            rng(seed_idx, "twister")
            x0 = randn(dimension, 1);
            if isfield(options, 'scale_x0')
                x0 = options.scale_x0 * x0;
            end
        end
        
        
        %% Test all selected methods for current seed
        for k_method = 1:length(method_names)
            mname = method_names{k_method};
            method_handle = selected_methods.(mname); % e.g., @SRD

            % Changed from params to options
            if strcmp(mname, 'BFGS')
                method_params = options.hanso;
                method_params.maxit = method_params.max_oracle_call;
                method_params.maxit_gradsamp = 0;
            elseif strcmp(mname, 'GradSamp')
                method_params = options.hanso;
                method_params.maxit = 0; 
                method_params.maxit_gradsamp = method_params.max_oracle_call; 
            else
                method_params = options.(mname);
            end
            
            fprintf('\n%s\n', repmat('-', 1, 50));
            fprintf('>>> Running Method: %s <<<\n', upper(mname));
            fprintf('%s\n', repmat('-', 1, 50));

            %% Setup logging
            log_filepath = local_setup_diary_logging(options, example, mname, dimension, seed_idx, x0_is_fixed);
            if ~isempty(log_filepath)
                diary(log_filepath); diary on;
                fprintf('Initial x0 (seed %d, method %s):\n', seed_idx, mname); disp(x0(1:min(10,end))');
            end
            
            start_time = tic;
            stats_this_run = struct();
            
            %% Run method
            try
                if contains(mname, 'SRD')
                    stats_this_run = method_handle(x0, problem_handles.func_value, problem_handles.G_oracle, problem_details.flags, method_params);
                elseif contains(mname, 'PBMDC')
                    stats_this_run = method_handle(x0, problem_handles.func_value1, problem_handles.func_value2, problem_handles.subgrad1, problem_handles.subgrad2, problem_details.flags, method_params);
                elseif strcmp(mname, 'BFGS') || strcmp(mname, 'GradSamp')
                    stats_this_run = hanso(x0, problem_handles.func_value, problem_handles.subgrad, problem_details.flags, method_params);
                else    % Run Polyak/NTD methods via Subgradient oracle using "run_optimizer"
                    if strcmp(mname, 'Polyak')
                        optimizer = Polyak(x0, problem_details.flags.opt_val);
                    elseif strcmp(mname, 'NTD')
                        optimizer = NTD(x0, inf, method_params);
                    else
                        error('Optimizer type "%s" does not have a defined execution path.', mname);
                    end
                    
                    [~, stats_this_run] = run_optimizer(optimizer, problem_handles.func_value, problem_handles.subgrad, problem_details.flags, method_params);
                end
                
                %% Collect results
                total_time = toc(start_time);
                stats_this_run.total_time = total_time;
                stats_this_run.seed_idx_ran = seed_idx;
                if strcmp(mname, 'SRD')
                    stats_this_run.name_of_method = 'SRDescent';
                elseif strcmp(mname, 'SRD_adapt')
                    stats_this_run.name_of_method = 'SRDescent-adapt';
                elseif strcmp(mname, 'NTD')
                    stats_this_run.name_of_method = 'NTDescent';
                else
                    stats_this_run.name_of_method = mname;
                end

                if isfield(options, 'outputs') && isfield(options.outputs, 'estimete_opt_val') && options.outputs.estimete_opt_val
                    stats_this_run.flags.opt_val = stats_this_run.objective_values_list(end);
                end
                
                % Aggregate results for this successful run
                aggregated_results.(mname).obj_vals(end+1) = stats_this_run.objective_values_list(end);
                aggregated_results.(mname).times(end+1) = total_time;
                aggregated_results.(mname).oracle_calls(end+1) = stats_this_run.num_oracle_list(end);
                aggregated_results.(mname).success_count = aggregated_results.(mname).success_count + 1;

                fprintf('\n%s: SUCCESS. Final Obj: %.4e, Time: %.2fs, Oracles: %d\n', upper(mname), ...
                        stats_this_run.objective_values_list(end), total_time, stats_this_run.num_oracle_list(end));


            catch ME 
                warning('\n%s: FAILURE for seed %d (%s)', mname, seed_idx, ME.message);
                aggregated_results.(mname).obj_vals(end+1) = NaN;
                aggregated_results.(mname).times(end+1) = NaN;
                aggregated_results.(mname).oracle_calls(end+1) = NaN;
                aggregated_results.(mname).failure_count = aggregated_results.(mname).failure_count + 1;
            end
            
            all_individual_run_stats{end+1} = stats_this_run;

            if ~isempty(log_filepath), diary off; end % Close diary for this specific run
        end
    end
    
    %% 6. Save Aggregated Numerical Results to .mat File
    local_save_mat_results(aggregated_results, example, options);


    %% 7. Generate Final Outputs (Plots, Tables) based on `options.outputs`
    % Inputs:
    %   statistics_list 
    %   log_plot       : (optional, default false) Boolean. If true, y-axis of some plots is logarithmic.
    %   contour_plot   : (optional) String specifying the contour plot type. For 2D,
    %                    use 'nesterov2D' or 'counter_example2D'; for 3D, 'nesterov3D-active-manifold'.
    % plot_oracle_vs_objective(example, statistics_list, true, [], false);
    results_summary = struct();
    results_summary.experiment_name = example.name;
    results_summary.aggregated_data = aggregated_results;
    
    output_base_dir = fullfile('results', options.experiment_group_name); % Use group name for subfolder

    if isfield(options, 'outputs') && isstruct(options.outputs)
        if isfield(options.outputs, 'generate_plots') && options.outputs.generate_plots
            if exist('plot_oracle_vs_objective', 'file')
                plot_oracle_vs_objective(all_individual_run_stats, ...
                    options.outputs.plot_log_y_axis, ... 
                    options.outputs.plot_contour_type);
                results_summary.plot_generation_attempted = true;
            else
                warning('Visualization function "plot_oracle_vs_objective.m" not found. Skipping plot generation.');
            end
        end

        if isfield(options.outputs, 'generate_table') && options.outputs.generate_table
            if exist('generate_table', 'file')
                table_output_dir = fullfile(output_base_dir, 'tables_txt');
                if ~exist(table_output_dir, 'dir'), mkdir(table_output_dir); end
                table_string = generate_table(example, aggregated_results, num_initials, options);
                results_summary.table_string = table_string;
                % Example: Save to file
                table_filename_txt = fullfile(table_output_dir, sprintf('%s_summary.txt', example.name));
                fid = -1; % Initialize file ID
                try
                    fid = fopen(table_filename_txt, 'w');
                    if fid == -1
                        error('generate_table:fopen_failed', 'Could not open file %s for writing. Check permissions or path.', table_filename_txt);
                    end
                    fprintf(fid, '%s', table_string);
                    fclose(fid);
                    fid = -1;
                    
                    results_summary.text_table_filepath = table_filename_txt;
                catch err_save_txt
                    if fid ~= -1
                        fclose(fid);
                    end
                    warning('Could not save readable text table to file: %s\nError ID: %s\nError Message: %s', ...
                            table_filename_txt, err_save_txt.identifier, err_save_txt.message);
                end
            end
        end
    end

    fprintf('\n%s\n', repmat('=', 1, 70));
    fprintf('===== Experiment Set Finished: %s =====\n', example.name);
    fprintf('%s\n\n', repmat('=', 1, 70));
end


% ----- Local Helper Functions for run_single_experiment_set.m -----

function [handles, details] = local_get_problem_handles_and_details(problem_func_name, example_struct, options_struct)
    handles = struct('func_value', [], 'subgrad', [], 'G_oracle', [], 'func_value1', [], 'func_value2', [], 'subgrad1', [], 'subgrad2', []);
    details = struct('flags', struct(), 'dimension', []); % Default dimension

    % Assumes problem functions are like:
    % [f, g, G, flags, d, f1, f2, g1, g2] = problem_name(example_struct_or_dimension_etc)
    % Or more simply: [f, g, G, flags, d] = problem_name(example_struct)
    problem_def_fh = str2func(problem_func_name);

    try
        if strcmp(problem_func_name, 'nesterov_ns2')
            
            [handles.func_value, handles.subgrad, handles.G_oracle, details.flags, ...
             handles.func_value1, handles.func_value2, handles.subgrad1, handles.subgrad2] = problem_def_fh(example_struct.dimension);
             details.dimension = example_struct.dimension;

        elseif strcmp(problem_func_name, 'max_of_smooth_degenerate_random')
            [handles.func_value, handles.subgrad, handles.G_oracle, details.flags, ...
             handles.func_value1, handles.func_value2, handles.subgrad1, handles.subgrad2] = ...
               problem_def_fh(example_struct.nb_functions, example_struct.dimension, example_struct.seed);
             details.dimension = example_struct.dimension;

        elseif strcmp(problem_func_name, 'min_of_smooth_random')
            [handles.func_value, handles.subgrad, handles.G_oracle, details.flags] = ...
                problem_def_fh(example_struct.nb_functions, example_struct.nb_features, example_struct.dimension, example_struct.seed);
            details.dimension = example_struct.dimension;

        elseif strcmp(problem_func_name, 'optval_LICQ')
             [handles.func_value, handles.subgrad, handles.G_oracle, details.flags] = ...
                problem_def_fh(example_struct.dimension_x, example_struct.dimension_y, example_struct.dimension_y_cons, example_struct.seed);
             details.dimension = example_struct.dimension_x;

        elseif strcmp(problem_func_name, 'counter_example_2D')
             [handles.func_value, handles.subgrad, handles.G_oracle, details.flags] = problem_def_fh();
             details.dimension = 2;
        end
    catch ME_prob_def
        error('Error calling problem definition function "%s": %s\nEnsure it is in src/problem_definitions/ and matches expected signature.', problem_func_name, ME_prob_def.message);
    end

    if ~isfield(details.flags, 'opt_val')
        if isfield(options_struct, 'default_opt_val')
            details.flags.opt_val = options_struct.default_opt_val;
        else
            details.flags.opt_val = Inf;
            % fprintf('Warning: Optimal value (flags.opt_val) not provided for problem %s. Using Inf.\n', problem_func_name);
        end
    end
end

function log_filepath = local_setup_diary_logging(options, example, mname, dimension, seed_idx, is_fixed_x0)
% Sets up diary logging for a single run if requested in options.
    log_filepath = '';
    if isfield(options, 'outputs') && isfield(options.outputs, 'enable_diary_log') && options.outputs.enable_diary_log
        log_base_dir = fullfile('results', options.experiment_group_name, 'logs', mname); % Group by method
        if ~exist(log_base_dir, 'dir')
            try mkdir(log_base_dir);
            catch err_mkdir
                warning('Could not create log directory: %s. Diary logging disabled for this run.\nError: %s', log_base_dir, err_mkdir.message);
                return;
            end
        end

        if is_fixed_x0
            log_filename = sprintf('%s_dim%d_%s_fixed_x0.txt', mname, dimension, example.name);
        else
            log_filename = sprintf('%s_dim%d_%s_seed%d.txt', mname, dimension, example.name, seed_idx);
        end
        log_filepath = fullfile(log_base_dir, log_filename);
    end
end

function local_save_mat_results(aggregated_results, example, options)
% Saves aggregated numerical results and individual stats to a .mat file.
    output_base_dir = fullfile('results', options.experiment_group_name);
    data_output_dir = fullfile(output_base_dir, 'data_matfiles');
    if ~exist(data_output_dir, 'dir'), mkdir(data_output_dir); end

    if isfield(options, 'x0_fixed') && ~isempty(options.x0_fixed)
        fn_suffix = '_fixed_x0';
    else
        fn_suffix = sprintf('_rand_seeds');
    end
    obj_tol_str = 'defaultTol';
    if isfield(options, 'SRD') && isfield(options.SRD, 'obj_tol') % Use a representative tolerance
        obj_tol_str = sprintf('%.0eTol', options.SRD.obj_tol);
    end

    save_filename = sprintf('results_%s_dim%d_%s%s.mat', example.name, aggregated_results.dimension, obj_tol_str, fn_suffix);
    full_save_path = fullfile(data_output_dir, save_filename);

    try
        save(full_save_path, 'aggregated_results', 'example', 'options');
    catch ME_save
        warning('Failed to save MAT results to %s: %s', full_save_path, ME_save.message);
    end
end
