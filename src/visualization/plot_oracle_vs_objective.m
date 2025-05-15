function figures = plot_oracle_vs_objective(statistics_list, log_plot, contour_plot)
%
% Inputs:
%   statistics_list : cell array of structs, one per method. Each struct should have:
%       - objective_values_list: [1 x N] numeric vector of f(x^k) values
%       - num_oracle_list      : [1 x N] numeric vector of cumulative oracle calls
%       - elapsed_time_list    : [1 x N] numeric vector of elapsed time values (optional)
%       - history              : [d x N] matrix of iterates (d = dimension)
%       - flags                : struct containing additional info (e.g., opt_val, opt_soln, params_for_legend)
%       - name_of_method       : string with the method name (for labeling)
%
%   log_plot       : (optional, default false) Boolean. If true, y-axis of some plots is logarithmic.
%   contour_plot   : (optional) String specifying the contour plot type. For 2D,
%                    use 'nesterov2D' or 'counter_example2D'; for 3D, 'nesterov3D-active-manifold'.
%
% Output:
%   figures : cell array containing figure and axis handles for each plot.
%
% Example:
%   figs = plot_oracle_vs_objective(statistics_list, true, 'nesterov2D', false);

% Set default values for optional arguments
if nargin < 2 || isempty(log_plot)
    log_plot = false;
end
if nargin < 3
    contour_plot = [];
end

% Set plotting parameters
markersize = 10;
linewidth = 2;
plot_marker_frequency = 1000; % (Not directly used; you may subsample if desired.)
scale_factor = 1; % Scale factor for quiver arrows in 3D

% Set default interpreter and font properties for all figures
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');
set(groot, 'DefaultAxesFontName', 'Times New Roman');
set(groot, 'DefaultAxesFontSize', 14);

%% Create base figures and axes
fig0 = figure('Position',[100 100 500 500]); % f(x^k) vs Oracle Calls
ax0 = axes(fig0);
fig1 = figure('Position',[100 100 500 500]); % f(x^k) vs Elapsed Time
ax1 = axes(fig1);

% Determine whether to show the distance plot (if flags.opt_soln exists)
last_flags = statistics_list{end}.flags;

% Determine the problem dimension from the last statistics entry.
if isfield(statistics_list{end}, 'x_result')
    dimension = length(statistics_list{end}.x_result);
else
    dimension = size(statistics_list{end}.history, 1);
end

if dimension == 2
    fig2 = figure('Position',[100 100 500 500]);
    ax2 = axes('Parent', fig2);
    global_x1_min = inf; global_x1_max = -inf;
    global_x2_min = inf; global_x2_max = -inf;
elseif dimension == 3
    fig2 = figure('Position',[100 100 600 600]);
    ax2 = axes('Parent', fig2, 'Projection', 'perspective');  % Use a valid projection
    global_x1_min = inf; global_x1_max = -inf;
    global_x2_min = inf; global_x2_max = -inf;
    global_x3_min = inf; global_x3_max = -inf;
else
    fig2 = [];
    ax2 = [];
end

% Get optimum value from flags (default to 0 if not present)
if isfield(last_flags, 'opt_val')
    opt_val = last_flags.opt_val;
else
    opt_val = 0;
end

%% Gather method names and update global ranges from histories
name_of_methods = {};
params_for_legends = {};
for i = 1:length(statistics_list)
    stats = statistics_list{i};
    name_of_methods{end+1} = stats.name_of_method;
    if isfield(stats.flags, 'params_for_legend')
        params_for_legends{end+1} = stats.flags.params_for_legend;
    else
        params_for_legends{end+1} = '';
    end
    history = stats.history;
    if dimension == 2
        global_x1_min = min(global_x1_min, min(history(1, :)));
        global_x1_max = max(global_x1_max, max(history(1, :)));
        global_x2_min = min(global_x2_min, min(history(2, :)));
        global_x2_max = max(global_x2_max, max(history(2, :)));
    elseif dimension == 3
        global_x1_min = min(global_x1_min, min(history(1, :)));
        global_x1_max = max(global_x1_max, max(history(1, :)));
        global_x2_min = min(global_x2_min, min(history(2, :)));
        global_x2_max = max(global_x2_max, max(history(2, :)));
        global_x3_min = min(global_x3_min, min(history(3, :)));
        global_x3_max = max(global_x3_max, max(history(3, :)));
    end
end

%% Create contour or surface plot (for 2D or 3D) if requested
if dimension == 2
    if ischar(contour_plot) && strcmp(contour_plot, 'nesterov2D')
        x1_range = global_x1_max - global_x1_min;
        x2_range = global_x2_max - global_x2_min;
        x_range = max(x1_range, x2_range);
        global_x1_min = global_x1_min - 1*x_range; %0.2
        global_x1_max = global_x1_max + 1*x_range;
        global_x2_min = global_x2_min - 1*x_range;
        global_x2_max = global_x2_max + 1*x_range;
        
        % Construct data points for contour
        x1_vals = linspace(global_x1_min, global_x1_max, 200);
        x2_vals = linspace(global_x2_min, global_x2_max, 200);
        [X1, X2] = meshgrid(x1_vals, x2_vals);
        F_vals = (1/4) * (X1 - 1).^2 + abs(X2 - 2*X1.^2 + 1);
        
        % Plot active manifold
        x1 = x1_vals;
        x2 = 2 * x1.^2 - 1; 
        hold(ax2, 'on');
        h = plot(ax2, x1, x2, 'Color', [1, 0.5, 0], 'LineWidth', 4, 'DisplayName', 'Active Manifold'); 
        set(h, 'Color', [1, 0.5, 0, 0.7]);  % Adding transparency using RGBA

        
        hold(ax2, 'on');
        % plot(ax2, 1, 1, 'p', 'MarkerSize', 15, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'DisplayName','Minimizer');
    elseif ischar(contour_plot) && strcmp(contour_plot, 'counter_example2D')
        x1_vals = linspace(global_x1_min-5, global_x1_max+5, 200);
        x2_vals = linspace(global_x2_min-10, global_x2_max+20, 200);
        [X1, X2] = meshgrid(x1_vals, x2_vals);
        F1 = -100 * ones(size(X1));
        F2 = 2*X1 + 3*X2;
        F3 = -2*X1 + 3*X2;
        F4 = 5*X1 + 2*X2;
        F5 = -5*X1 + 2*X2;
        F_vals = max(cat(3, F1, F2, F3, F4, F5), [], 3);
        
        hold(ax2, 'on');
        % plot(ax2, 0, -50, 'p', 'MarkerSize', 15, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'DisplayName', 'Minimizer');
    elseif ~isempty(contour_plot)
        error("Unknown contour_plot value. Please use 'nesterov2D' or 'counter_example2D'.");
    end
    
    % Plot contour
    if ~isempty(contour_plot) 
        hold(ax2, 'on');
        contourf(ax2, X1, X2, F_vals, 10, 'LineColor','none', 'FaceAlpha', 0.3, 'HandleVisibility', 'off');
        contour(ax2, X1, X2, F_vals, 10, 'LineWidth', 0.8, 'FaceAlpha', 0.3, 'LineStyle','-', 'LineColor','k', 'HandleVisibility', 'off');
    end
elseif dimension == 3
    if ischar(contour_plot) && strcmp(contour_plot, 'nesterov3D-active-manifold')
        x1_range = global_x1_max - global_x1_min;
        global_x1_min = global_x1_min - 0.2*x1_range;
        global_x1_max = global_x1_max + 0.2*x1_range;
        x1 = linspace(global_x1_min, global_x1_max, 300);
        x2 = 2 * x1.^2 - 1; 
        x3 = 2 * x2.^2 - 1; 
        hold(ax2, 'on');
        plot3(ax2, x1, x2, x3, 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'Active Manifold'); 
        % surf(ax2, X1, X2, X3, 'FaceColor', [0.678, 0.847, 0.902], 'LineWidth', 1.5, 'FaceAlpha', 0.5, 'EdgeColor','k', 'DisplayName', 'Active Manifold');
        hold(ax2, 'on');
        plot3(ax2, 1, 1, 1, 'p', 'MarkerSize', 15, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'DisplayName', 'Minimizer');
    elseif ~isempty(contour_plot)
        error("Unknown contour_plot value. For 3D use 'nesterov3D-active-manifold'.");
    end
end


%% Set up color and line style mapping.
unique_methods = unique(name_of_methods);
colors = lines(numel(unique_methods)); % Use MATLAB's lines colormap
color_mapping = containers.Map();
for i = 1:numel(unique_methods)
    % color_mapping(unique_methods{i}) = colors(i,:);
    if strcmp(unique_methods{i}, 'Polyak')
        color_mapping(unique_methods{i}) = [0 0 0];  % solid line
    elseif strcmp(unique_methods{i}, 'SRDescent')
        color_mapping(unique_methods{i}) = [0.9290 0.6940 0.1250];  % solid line
    elseif strcmp(unique_methods{i}, 'SRDescent-adapt')
        color_mapping(unique_methods{i}) = [0.8500 0.3250 0.0980]; % dashed line
    elseif strcmp(unique_methods{i}, 'Bundle')
        color_mapping(unique_methods{i}) = [0.4940 0.1840 0.5560]; % dash-dotted line
    elseif strcmp(unique_methods{i}, 'GradSamp')
        color_mapping(unique_methods{i}) = [0 0.4470 0.7410];
    elseif strcmp(unique_methods{i}, 'BFGS')
        color_mapping(unique_methods{i}) = [0.5 0.5 0.5];
    elseif strcmp(unique_methods{i}, 'NTDescent')
        color_mapping(unique_methods{i}) = [0.3010 0.7450 0.9330];
    else
        color_mapping(unique_methods{i}) = colors(i,:);
    end
end

% Define line style (dashes). In MATLAB, use: '-' (solid), '--' (dashed), ':' (dotted)
dashes_dict = containers.Map();
for i = 1:numel(unique_methods)
    if strcmp(unique_methods{i}, 'SRDescent') || strcmp(unique_methods{i}, 'Polyak')
        dashes_dict(unique_methods{i}) = '-';  % solid line
    elseif strcmp(unique_methods{i}, 'SRDescent-adapt')
        dashes_dict(unique_methods{i}) = '--'; % dashed line
    elseif strcmp(unique_methods{i}, 'Bundle') || strcmp(unique_methods{i}, 'NTDescent')
        dashes_dict(unique_methods{i}) = '-.'; % dash-dotted line
    elseif strcmp(unique_methods{i}, 'BFGS')
        dashes_dict(unique_methods{i}) = ':';  % dotted line
    else
        dashes_dict(unique_methods{i}) = '--';  % dotted line
    end
end

%% Plot curves for each method.
for i = 1:length(statistics_list)
    stats = statistics_list{i};
    obj_vals = stats.objective_values_list;
    oracle_calls = stats.num_oracle_list;
    if isfield(stats, 'elapsed_time_list')
        elapsed_time = stats.elapsed_time_list;
    else
        elapsed_time = [];
    end
    history = stats.history;
    method_name = stats.name_of_method;
    if isfield(stats.flags, 'params_for_legend')
        params_for_legend = stats.flags.params_for_legend;
    else
        params_for_legend = '';
    end
    
    % Process objective values for certain methods if needed.
    if strcmp(method_name, 'PolyakSGM')
        obj_vals = cummin(obj_vals);
    end
    
    % f(x^k) - opt_val vs Oracle Calls on ax0.
    y_vals = obj_vals - opt_val;
    plot(ax0, oracle_calls, y_vals, 'LineWidth', linewidth, ...
        'Color', color_mapping(method_name), 'LineStyle', dashes_dict(method_name), ...
        'DisplayName', method_name);
    hold(ax0, 'on');
    
    % f(x^k) - opt_val vs Elapsed Time on ax1 (if available).
    if ~isempty(elapsed_time)
        plot(ax1, elapsed_time, y_vals, 'LineWidth', linewidth, ...
            'Color', color_mapping(method_name), 'LineStyle', dashes_dict(method_name), ...
            'DisplayName', method_name);
        hold(ax1, 'on');
    end
    
    % Trajectory plot on ax2.
    if dimension == 2 && ~isempty(ax2)
        X = history(1,1:end-1);
        Y = history(2,1:end-1);
        U = diff(history(1,:));
        V = diff(history(2,:)); 
        quiver(ax2, X, Y, U, V, 0, 'Color', color_mapping(method_name), 'LineWidth', linewidth, ...
            'DisplayName', method_name, 'MaxHeadSize', 0.15);
        hold(ax2, 'on');
        scatter(ax2, history(1,:), history(2,:), 40, color_mapping(method_name), 'filled', 'HandleVisibility', 'off');
        % scatter(ax2, history(1,:), history(2,:), 50, 'k', 'filled', 'DisplayName', 'Iterates of Gradient Sampling');
    elseif dimension == 3 && ~isempty(ax2)
        X = history(1,1:end-1);
        Y = history(2,1:end-1);
        Z = history(3,1:end-1);
        U = scale_factor * diff(history(1, :));
        V = scale_factor * diff(history(2, :));
        W = scale_factor * diff(history(3, :)); 
        quiver3(ax2, X, Y, Z, U, V, W, 0, 'Color', color_mapping(method_name), 'LineWidth', linewidth, ...
            'DisplayName', method_name, 'MaxHeadSize', 0.15);
        hold(ax2, 'on');
        scatter3(ax2, history(1,:), history(2,:), history(3,:), 17, color_mapping(method_name), 'filled', 'HandleVisibility', 'off');
    end
end

if log_plot
    set(ax0, 'YScale', 'log');
    set(ax1, 'YScale', 'log');
end

% Set axis labels and legends.
xlabel(ax0, 'Cumulative Oracle Calls');
ylabel(ax0, '$f(x^k)-f^\ast$', 'Interpreter', 'latex');
legend(ax0, 'Location', 'northeast');

xlabel(ax1, 'Elapsed Time (s)');
ylabel(ax1, '$f(x^k)-f^\ast$', 'Interpreter', 'latex');
legend(ax1, 'Location', 'northeast');

if dimension == 2 && ~isempty(ax2)
    xlabel(ax2, '$x_1$', 'Interpreter', 'latex');
    ylabel(ax2, '$x_2$', 'Interpreter', 'latex');
    title(ax2, 'Trajectory of Iterates', 'Interpreter', 'latex');
    legend(ax2, 'show');
elseif dimension == 3 && ~isempty(ax2)
    % Add grid
    grid on;
    xlabel(ax2, '$x_1$', 'Interpreter', 'latex');
    ylabel(ax2, '$x_2$', 'Interpreter', 'latex');
    zlabel(ax2, '$x_3$', 'Interpreter', 'latex');
    title(ax2, 'Trajectory of Iterates', 'Interpreter', 'latex');
    legend(ax2, 'show');
end

drawnow;

% Collect all figure and axis handles into a cell array for return.
figures = {fig0, ax0, fig1, ax1};
if ~isempty(ax2)
    figures = [figures, {fig2, ax2}];
end

end
