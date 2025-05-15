function table_string = generate_table(example, results, num_initials, options)
    output_str = '';

    if isfield(options, 'methods_to_run') && ~isempty(options.methods_to_run)
        for k_idx = 1:length(options.methods_to_run)
            method_name_temp = options.methods_to_run{k_idx};
            if strcmpi(method_name_temp, 'BFGS') || strcmpi(method_name_temp, 'GradSamp')
                if isfield(options, 'hanso') && isfield(options.hanso, 'obj_tol')
                    obj_tol = options.hanso.obj_tol;
                    break;
                end
            % Check if the method_name_temp itself is a field in options
            elseif isfield(options, method_name_temp) 
                if isfield(options.(method_name_temp), 'obj_tol')
                    obj_tol = options.(method_name_temp).obj_tol;
                    break;
                end
            end
        end
    end

    % Define column widths
    dim_col_w = 6;        % For "Dim." and its values
    method_col_w = 18;    % For "Method" and method names
    failures_col_w = 10;  % For "Failures"
    time_col_w = 22;      % For "Time (s)" and "mean (failed_mean)"
    obj_val_col_w = 22;   % For "Obj Val"
    oracles_col_w = 22;   % For "Oracles"

    % Table Title and Context
    output_str = [output_str sprintf('\nComparison on %s with %d initial points.\n', example.name, num_initials)];

    header_str = sprintf('| %-*s | %-*s | %*s | %-*s | %-*s | %-*s |', ...
        dim_col_w, 'Dim.', ...
        method_col_w, 'Method', ...
        failures_col_w, 'Failures', ...
        time_col_w, 'Time (s)', ...
        obj_val_col_w, 'Obj Val', ...
        oracles_col_w, 'Oracles');
    
    separator_line = repmat('-', 1, length(header_str));
    
    output_str = [output_str sprintf('%s\n', separator_line)];
    output_str = [output_str sprintf('%s\n', header_str)];
    output_str = [output_str sprintf('%s\n', separator_line)];
    
    if ~isfield(results, 'dimension')
        padded_message = pad('Results dimension not provided.', length(header_str)-2, 'center');
        output_str = [output_str sprintf('| %s |\n', padded_message)];
        output_str = [output_str sprintf('%s\n', separator_line)];
        table_string = output_str;
        return;
    end
    dim = results.dimension;
    
    if ~isfield(options, 'methods_to_run') || isempty(options.methods_to_run)
        padded_message = pad('No methods specified in options.methods_to_run.', length(header_str)-2, 'center');
        output_str = [output_str sprintf('| %s |\n', padded_message)];
        output_str = [output_str sprintf('%s\n', separator_line)];
        table_string = output_str;
        return;
    end
    method_names = options.methods_to_run;
    
    % Format string for error rows (uses dynamic width specifiers '*')
    error_message_span_width = failures_col_w + time_col_w + obj_val_col_w + oracles_col_w + 3*3; % 3 separators " | "
    error_row_format_spec = '| %*d | %-*s | %-*s |\n';

    for m = 1:length(method_names)
        mname = method_names{m};
        
        if ~isfield(results, mname) || ~isstruct(results.(mname)) || ...
           ~isfield(results.(mname), 'obj_vals') || ~isfield(results.(mname), 'times') || ~isfield(results.(mname), 'oracle_calls')
            error_text = sprintf('Data missing or incomplete for method: %s.', mname);
            padded_error_text = pad(error_text, error_message_span_width, 'center');
            output_str = [output_str sprintf(error_row_format_spec, ...
                            dim_col_w, dim, ...
                            method_col_w, mname, ...
                            error_message_span_width, padded_error_text)];
            continue; % Skip to next method
        end
        data = results.(mname);
        
        current_method_obj_tol = NaN;
        if strcmpi(mname, 'BFGS') || strcmpi(mname, 'GradSamp')
            if isfield(options, 'hanso') && isfield(options.hanso, 'obj_tol')
                current_method_obj_tol = options.hanso.obj_tol;
            end
        elseif isfield(options, mname) && isfield(options.(mname), 'obj_tol') % Check options.mname.obj_tol
             current_method_obj_tol = options.(mname).obj_tol;
        end
        
        if isnan(current_method_obj_tol)
            error_text = sprintf('Obj_tol not defined for method %s. Cannot compute stats.', mname);
            padded_error_text = pad(error_text, error_message_span_width, 'center');
            output_str = [output_str sprintf(error_row_format_spec, ...
                            dim_col_w, dim, ...
                            method_col_w, mname, ...
                            error_message_span_width, padded_error_text)];
            continue;
        end
        
        summary = calculate_stats(...
            data.obj_vals, data.times, data.oracle_calls, ...
            current_method_obj_tol, num_initials);
        
        data_row_str = sprintf('| %*d | %-*s | %*d | %-*s | %-*s | %-*s |\n', ...
            dim_col_w, dim, ...
            method_col_w, mname, ...
            failures_col_w, summary.failures, ...
            time_col_w, summary.time, ...
            obj_val_col_w, summary.obj_val, ...
            oracles_col_w, summary.oracles);
        output_str = [output_str data_row_str];
    end
    
    output_str = [output_str sprintf('%s\n', separator_line)];
    output_str = [output_str newline];
    
    table_string = output_str;
end

% --- Helper Functions (from the original script) ---

function summary = calculate_stats(objs, times, oracles, obj_tol, n)
%calculate_stats Compute statistics
    success = (objs <= obj_tol); % Ensure objs and obj_tol are comparable
    summary.failures = n - sum(success);
    
    % Handle cases where all success or all fail for times(success) etc.
    summary.time = format_mean(times(success), '%.1e');
    summary.obj_val = format_mean(objs(success), '%.1e');
    summary.oracles = format_mean(oracles(success), '%.1e');
    
    if summary.failures > 0 
        % Only add stats for non-successful runs if there are any
        if any(~success)
            summary.time = sprintf('%s (%s)', summary.time, format_mean(times(~success), '%.1e')); 
            summary.obj_val = sprintf('%s (%s)', summary.obj_val, format_mean(objs(~success), '%.1e'));
            summary.oracles = sprintf('%s (%s)', summary.oracles, format_mean(oracles(~success), '%.1e'));
        end
    end
end

function str = format_mean(values, fmt)
    if ~isnumeric(values)
        str = '--';
        warning('format_mean: Input values are not numeric.');
        return;
    end

    valid_values = values(~isnan(values)); 
    
    if isempty(valid_values)
        str = '--';
    else
        str = sprintf(fmt, mean(valid_values));
        str = regexprep(str, 'e([+-])0(\d)$', 'e$1$2');
    end
end