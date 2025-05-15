clear; clc; close all;

experiment_group_name = 'group3_min_of_smooth';

%% Define example parameters
example = struct();
example.name         = 'min_of_smooth_random';
example.dimension    = 300;
example.nb_functions = 10;  % 50 100
example.nb_features  = 300;
example.seed = 123;      % seed for generating random instances

%% Create options with custom parameters
options = struct();
options.experiment_group_name = experiment_group_name; % For results folder

% Specify which methods to run
options.methods_to_run = {'SRD', 'SRD_adapt', 'Polyak', 'GradSamp', 'NTD', 'BFGS'};

%%%%%%%%% Control outputs %%%%%%%%%
options.outputs.enable_diary_log  = true;      % Enable/disable diary logs
options.outputs.generate_plots    = true;      % Generate plots using plot_oracle_vs_objective
options.outputs.estimete_opt_val  = false;     % Estimate optimal value using the final objective
options.outputs.plot_log_y_axis   = true;      % Custom flag for your plotting function
options.outputs.plot_contour_type = [];        % 'nesterov2D' or 'counter_example2D' for 2D; 'nesterov3D-active-manifold' for 3D
options.outputs.generate_table    = false;      % Generate table summary

%%%%%%%%% Printing & Log Info %%%%%%%%%
options.Polyak.print_frequency = 100000;       % per oracle
options.PBMDC.print_frequency  = 10;           % per serious step
options.NTD.print_frequency    = 100;          % per outer iter
options.SRD.print_frequency    = 100;          % per outer iter
options.hanso.print_frequency  = 10;           % per outer iter

options.Polyak.save_frequency = 1; 
options.NTD.save_frequency = 1; 

%%%%%%%%% Termination Criteria %%%%%%%%%
options.scale_x0 = 1;

% options.SRD.alpha = 0.5; (default)
% options.SRD.beta = 1e-4; (default)

options.Polyak.max_oracle_call = Inf;
options.PBMDC.max_oracle_call = Inf;
options.NTD.max_oracle_call = Inf;
options.SRD.max_oracle_call = Inf;
options.hanso.max_oracle_call = Inf;

options.Polyak.max_time_seconds = 60;
options.PBMDC.max_time_seconds = 60;
options.NTD.max_time_seconds = 60;
options.SRD.max_time_seconds = 60;
options.hanso.max_time_seconds = 60;

options.SRD.nu_tol = 0;
options.SRD.epsilon_tol = 0;

options.Polyak.obj_tol = 1e-10;
options.PBMDC.obj_tol = 1e-10;
options.NTD.obj_tol = 1e-10;
options.SRD.obj_tol = 1e-10;
options.hanso.obj_tol = 1e-10;

options.hanso.normtol = 0;          % stopping of BFGS
% options.hanso.evaldist = 1e-4;    % both for BFGS & GS
% options.hanso.eps_opt = 1e-4;     % GS termination tolerance (default: 1e-6)
% options.hanso.nu_opt = 1e-4;      % GS stationarity termination tolerance (default: 1e-6)
% options.hanso.delta = 1e-16;      % tolerance for line search (default: 1e-16)

options.PBMDC.step_tol = 0;


%% test for a fixed/random initial points
% x0 = (-1).^(1:example.dimension) * 0.5;
% options.x0 = x0(:); % Converts to a column vector

num_initials = 1;   % If options.x0_fixed is set, this will effectively be 1

results_summary = run_single_experiment(example, num_initials, options);


%% Post-processing or Display Summary
if isfield(results_summary, 'table_string') && ~isempty(results_summary.table_string)
    disp(results_summary.table_string);
    if isfield(results_summary, 'text_table_filepath')
        fprintf('INFO: Summary table saved to: %s\n', results_summary.text_table_filepath);
    end
end
