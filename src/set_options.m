function options = set_options(options)
%SET_OPTIONS Set and validate algorithm options with defaults
%   Validates input options and sets default values for missing parameters
%   for SRD, Polyak, and NTD algorithms. Throws errors for invalid inputs.

    if nargin < 1
        options = struct();
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%% Polyak Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isfield(options, 'Polyak')
        options.Polyak = struct();
    end
    
    % c1 parameter
    % if ~isfield(options.Polyak, 'c1')
    %     options.Polyak.c1 = 0;
    % else
    %     validate_nonnegative_real(options.Polyak.c1, 'Polyak.c1');
    % end
    
    % c2 parameter
    % if ~isfield(options.Polyak, 'c2')
    %     options.Polyak.c2 = 0.5;
    % else
    %     validate_scalar_in_interval(options.Polyak.c2, 'Polyak.c2', 0, 1);
    % end
    
    % Max oracle calls
    if ~isfield(options.Polyak, 'max_oracle_call')
        options.Polyak.max_oracle_call = 2e6;
    else
        validate_positive_integer(options.Polyak.max_oracle_call, 'Polyak.max_oracle_call');
    end
    
    % Time limit
    if ~isfield(options.Polyak, 'max_time_seconds')
        options.Polyak.max_time_seconds = 60;
    else
        validate_positive_real(options.Polyak.max_time_seconds, 'Polyak.max_time_seconds');
    end
    
    % Print frequency
    if ~isfield(options.Polyak, 'print_frequency')
        options.Polyak.print_frequency = 100000;
    else
        validate_positive_integer(options.Polyak.print_frequency, 'Polyak.print_frequency');
    end
    
    % Save frequency (prevent out of memory due to large iter num)
    if ~isfield(options.Polyak, 'save_frequency')
        options.Polyak.print_frequency = 1000;
    else
        validate_positive_integer(options.Polyak.save_frequency, 'Polyak.save_frequency');
    end

    % Objective tolerance
    if ~isfield(options.Polyak, 'obj_tol')
        options.Polyak.obj_tol = 1e-6;
    else
        validate_positive_real(options.Polyak.obj_tol, 'Polyak.obj_tol');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% NTD Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isfield(options, 'NTD')
        options.NTD = struct();
    end
    
    % Adaptive grid size
    if ~isfield(options.NTD, 'adaptive_grid_size')
        options.NTD.adaptive_grid_size = false;
    else
        validate_logical_scalar(options.NTD.adaptive_grid_size, 'NTD.adaptive_grid_size');
    end
    
    % Trust region usage
    if ~isfield(options.NTD, 'use_trust_region')
        options.NTD.use_trust_region = true;
    else
        validate_logical_scalar(options.NTD.use_trust_region, 'NTD.use_trust_region');
    end
    
    % Scale factor
    if ~isfield(options.NTD, 's_scale_factor')
        options.NTD.s_scale_factor = 1e-6;
    else
        validate_positive_real(options.NTD.s_scale_factor, 'NTD.s_scale_factor');
    end
    
    % Max oracle calls
    if ~isfield(options.NTD, 'max_oracle_call')
        options.NTD.max_oracle_call = 2e6;
    else
        validate_positive_integer(options.NTD.max_oracle_call, 'NTD.max_oracle_call');
    end
    
    % Time limit
    if ~isfield(options.NTD, 'max_time_seconds')
        options.NTD.max_time_seconds = 60;
    else
        validate_positive_real(options.NTD.max_time_seconds, 'NTD.max_time_seconds');
    end
    
    % Print frequency
    if ~isfield(options.NTD, 'print_frequency')
        options.NTD.print_frequency = 100000;
    else
        validate_positive_integer(options.NTD.print_frequency, 'NTD.print_frequency');
    end
    
    % Save frequency (prevent out of memory due to large iter num)
    if ~isfield(options.NTD, 'save_frequency')
        options.NTD.print_frequency = 1;
    else
        validate_positive_integer(options.NTD.save_frequency, 'NTD.save_frequency');
    end

    % Objective tolerance
    if ~isfield(options.NTD, 'obj_tol')
        options.NTD.obj_tol = 1e-6;
    else
        validate_positive_real(options.NTD.obj_tol, 'NTD.obj_tol');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% SRD Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isfield(options, 'SRD')
        options.SRD = struct();
    end
    
    % Initial regularization parameter
    if ~isfield(options.SRD, 'epsilon_0_0')
        options.SRD.epsilon_0_0 = 5;
    else
        validate_positive_real(options.SRD.epsilon_0_0, 'SRD.epsilon_0_0');
    end

    % Epsilon tolerance
    if ~isfield(options.SRD, 'epsilon_tol')
        options.SRD.epsilon_tol = 1e-2;
    else
        validate_nonnegative_real(options.SRD.epsilon_tol, 'SRD.epsilon_tol');
    end

    % Initial stationary target
    if ~isfield(options.SRD, 'nu_0')
        options.SRD.nu_0 = 1e-2;
    else
        validate_positive_real(options.SRD.nu_0, 'SRD.nu_0');
    end

    % Stationary tolerance
    if ~isfield(options.SRD, 'nu_tol')
        options.SRD.nu_tol = 1e-4;
    else
        validate_nonnegative_real(options.SRD.nu_tol, 'SRD.nu_tol');
    end
    
    % Reduction factor for epsilon
    if ~isfield(options.SRD, 'alpha_epsilon')
        options.SRD.alpha_epsilon = 0.9;
    elseif options.SRD.alpha_epsilon >= 1 || options.SRD.alpha_epsilon <= 0
        error('set_options_SRD: Reduction factor alpha_epsilon must be in (0, 1)');
    end

    % Reduction factor for nu
    if ~isfield(options.SRD, 'alpha_nu')
        options.SRD.alpha_nu = 0.5;
    elseif options.SRD.alpha_nu >= 1 || options.SRD.alpha_nu <= 0
        error('set_options_SRD: Reduction factor alpha_nu must be in (0, 1)');
    end

    % Constant in Armijo's line search
    if ~isfield(options.SRD, 'beta')
        options.SRD.beta = 1e-4;
    elseif options.SRD.beta >= 1 || options.SRD.beta <= 0
        error('set_options_SRD: Armijo-like parameter beta must be in (0, 1)');
    end

    % Objective tolerance (optional)
    if ~isfield(options.SRD, 'obj_tol')
        options.SRD.obj_tol = 1e-6;
    else
        validate_nonnegative_real(options.SRD.obj_tol, 'SRD.obj_tol');
    end
    
    % Maximum oracle calls
    if ~isfield(options.SRD, 'max_oracle_call')
        options.SRD.max_oracle_call = 1e6;
    else
        validate_positive_integer(options.SRD.max_oracle_call, 'SRD.max_oracle_call');
    end
    
    % Maximum time
    if ~isfield(options.SRD, 'max_time_seconds')
        options.SRD.max_time_seconds = 60;
    else
        validate_positive_real(options.SRD.max_time_seconds, 'SRD.max_time_seconds');
    end
    
    % Print frequency
    if ~isfield(options.SRD, 'print_frequency')
        options.SRD.print_frequency = 1000;
    else
        validate_positive_integer(options.SRD.print_frequency, 'SRD.print_frequency');
    end
    
    % Save frequency 
    if ~isfield(options.SRD, 'save_frequency')
        options.SRD.save_frequency = 1;
    else
        validate_positive_integer(options.SRD.save_frequency, 'SRD.save_frequency');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%% SRD_adapt Options %%%%%%%%%%%%%%%%%%%%%%%%
    options.SRD_adapt = options.SRD;

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%% RedistBundle Options %%%%%%%%%%%%%%%%%%%%%%%%
    % % Convexification growth parameter
    % if ~isfield(options.RedistBundle, 'Gamma')
    %     options.RedistBundle.Gamma = 4;
    % elseif Gamma <= 1
    %     error('set_options: Gamma must be larger than 1');
    % end
    % 
    % % Unacceptable increase parameter
    % if ~isfield(options.RedistBundle, 'M0')
    %     options.RedistBundle.M0 = 10;
    % else
    %     validate_positive_real(options.RedistBundle.M0, 'RedistBundle.M0');
    % end
    % 
    % % Initial R parameter
    % if ~isfield(options.RedistBundle, 'R0')
    %     options.RedistBundle.R0 = 1;
    % else
    %     validate_positive_real(options.RedistBundle.R0, 'RedistBundle.R0');
    % end
    % 
    % % Armijo-like parameter
    % if ~isfield(options.RedistBundle, 'm')
    %     options.RedistBundle.m = 0.01;
    % elseif m >= 1 || m <= 0
    %     error('set_options: Armijo-like parameter m must be in (0, 1)');
    % end
    % 
    % % Selection strategy
    % % 1: I = [1,2,...,n+1]               (size = n+1)
    % % 2: I = J_active & [n+1, i]         (size < n+1)
    % % 3: I = [n+1,i] & [0(aggregated)]   (size =  3 )
    % if ~isfield(options.RedistBundle, 'strat')
    %     options.RedistBundle.strat = 3;
    % elseif ~ismember(options.RedistBundle.strat, [1, 2, 3])
    %     error('set_options: strat must be in [1, 2, 3]');
    % end
    % 
    % % Delta tolerance
    % if ~isfield(options.RedistBundle, 'delta_tol')
    %     options.RedistBundle.delta_tol = 1e-6;
    % else
    %     validate_nonnegative_real(options.RedistBundle.delta_tol, 'RedistBundle.delta_tol');
    % end
    % 
    % % Objective tolerance (optional)
    % if ~isfield(options.RedistBundle, 'obj_tol')
    %     options.RedistBundle.obj_tol = 1e-6;
    % else
    %     validate_nonnegative_real(options.RedistBundle.obj_tol, 'RedistBundle.obj_tol');
    % end
    % 
    % % Maximum oracle calls
    % if ~isfield(options.RedistBundle, 'max_oracle_call')
    %     options.RedistBundle.max_oracle_call = 1e6;
    % else
    %     validate_positive_integer(options.RedistBundle.max_oracle_call, 'RedistBundle.max_oracle_call');
    % end
    % 
    % % Maximum time
    % if ~isfield(options.RedistBundle, 'max_time_seconds')
    %     options.RedistBundle.max_time_seconds = 60;
    % else
    %     validate_positive_real(options.RedistBundle.max_time_seconds, 'RedistBundle.max_time_seconds');
    % end
    % 
    % % Print frequency
    % if ~isfield(options.RedistBundle, 'print_frequency')
    %     options.RedistBundle.print_frequency = 1000;
    % else
    %     validate_positive_integer(options.RedistBundle.print_frequency, 'RedistBundle.print_frequency');
    % end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% PBMDC Options %%%%%%%%%%%%%%%%%%%%%%%%
    if ~isfield(options.PBMDC, 'QPdual')
        options.PBMDC.QPdual = true;
    elseif ~ismember(options.PBMDC.QPdual, [true, false])
        error('set_options: PBMDC.QPdual must be true/false');
    end

    if ~isfield(options.PBMDC, 'kappa')
        options.PBMDC.kappa = 0.1;
    elseif options.PBMDC.kappa >= 1 || options.PBMDC.kappa <= 0
        error('set_options: PBMDC.kappa must be in (0,1)');
    end
    
    if ~isfield(options.PBMDC, 'tmin')
        options.PBMDC.tmin = 1e-5;
    else
        validate_positive_real(options.PBMDC.tmin, 'PBMDC.tmin');
    end

    if ~isfield(options.PBMDC, 't')
        options.PBMDC.t = 1;
    else
        validate_positive_real(options.PBMDC.t, 'PBMDC.t');
    end

    if ~isfield(options.PBMDC, 'tmax')
        options.PBMDC.tmax = 1e5;
    else
        validate_positive_real(options.PBMDC.tmax, 'PBMDC.tmax');
    end

    if ~isfield(options.PBMDC, 'nBun')
        options.PBMDC.nBun = 100; %max(100,min(n+5,1000));
    else
        validate_positive_integer(options.PBMDC.nBun, 'PBMDC.nBun');
    end

    if ~isfield(options.PBMDC, 'nBun2')
        options.PBMDC.nBun2 = 1;
    else
        validate_positive_integer(options.PBMDC.nBun2, 'PBMDC.nBun2');
    end

    % Solving the QP by Matlab quadprog (True) or Gurobi (False)
    if ~isfield(options.PBMDC, 'matlab')
        options.PBMDC.matlab = true;
    elseif ~ismember(options.PBMDC.matlab, [true, false])
        error('set_options: PBMDC.matlab must be true/false');
    end

    % Step tolerance (||trial point - prox center||)
    if ~isfield(options.PBMDC, 'step_tol')
        options.PBMDC.step_tol = 1e-4;
    else
        validate_nonnegative_real(options.PBMDC.step_tol, 'PBMDC.step_tol');
    end

    % Objective tolerance (optional)
    if ~isfield(options.PBMDC, 'obj_tol')
        options.PBMDC.obj_tol = 1e-6;
    else
        validate_nonnegative_real(options.PBMDC.obj_tol, 'PBMDC.obj_tol');
    end
    
    % Maximum (serious) iteration number
    if ~isfield(options.PBMDC, 'max_iter')
        options.PBMDC.max_iter = 5000;
    else
        validate_positive_integer(options.PBMDC.max_iter, 'PBMDC.max_iter');
    end

    % Maximum oracle calls
    if ~isfield(options.PBMDC, 'max_oracle_call')
        options.PBMDC.max_oracle_call = 1e6;
    else
        validate_positive_integer(options.PBMDC.max_oracle_call, 'PBMDC.max_oracle_call');
    end
    
    % Maximum time
    if ~isfield(options.PBMDC, 'max_time_seconds')
        options.PBMDC.max_time_seconds = 60;
    else
        validate_positive_real(options.PBMDC.max_time_seconds, 'PBMDC.max_time_seconds');
    end
    
    % Print frequency
    if ~isfield(options.PBMDC, 'print_frequency')
        options.PBMDC.print_frequency = 10;
    else
        validate_positive_integer(options.PBMDC.print_frequency, 'PBMDC.print_frequency');
    end

    %%%%%%%%%%%%%%%%%%%%%%% Deprecated Parameters Check %%%%%%%%%%%%%%%%%%%%
    check_deprecated(options, 'SRD', {});
    check_deprecated(options, 'Polyak', {});
    check_deprecated(options, 'NTD', {});
end

%% Validation Helper Functions
function validate_positive_real(value, name)
    if ~(isscalar(value) && isreal(value) && value > 0)
        error('set_options: %s must be a positive real scalar', name);
    end
end

function validate_nonnegative_real(value, name)
    if ~(isscalar(value) && isreal(value) && value >= 0)
        error('set_options: %s must be a non-negative real scalar', name);
    end
end

function validate_positive_integer(value, name)
    if ~(isscalar(value) || value <= 0 || rem(value,1) ~= 0)
        error('set_options: %s must be a positive integer', name);
    end
end

function validate_scalar_in_interval(value, name, lb, ub)
    if ~(isscalar(value) && isreal(value) && value > lb && value < ub)
        error('set_options: %s must be in interval (%g, %g)', name, lb, ub);
    end
end

function validate_logical_scalar(value, name)
    if ~(islogical(value) && isscalar(value))
        error('set_options: %s must be a logical scalar', name);
    end
end

function validate_nonnegative_real_or_empty(value, name)
    if ~isempty(value) && ~(isscalar(value) && isreal(value) && value >= 0)
        error('set_options: %s must be non-negative real or empty', name);
    end
end

function check_deprecated(options, field, deprecated_list)
    if isfield(options, field)
        fnames = fieldnames(options.(field));
        for i = 1:length(fnames)
            if any(strcmp(fnames{i}, deprecated_list))
                error('set_options: %s.%s is deprecated', field, fnames{i});
            end
        end
    end
end