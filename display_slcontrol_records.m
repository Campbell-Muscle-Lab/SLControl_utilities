function output=display_slcontrol_records(varargin)

params.record_file_strings={};
params.slr_pCa_values=[];
params.force_subplot=[];
params.fl_subplot=[];
params.sl_subplot=[];
params.start_time_s=[];
params.stop_time_s=[];
params.t_offset_s=0;
params.normalize_fl=1;
params.max_points_per_line=5000;
params.trace_line_width=1;
params.trace_colors=[1 0 0;0 1 0;0 0 1];
params.single_fl_color=[];
params.display_pCa_values=1;
params.pCa_font_size=12;
params.pCa_x_offset=1.15;
params.pCa_label_y_offset=1.2;
params.force_error_label=''; % options are 'sem' or 'std'
params.patch_color=0.75*ones(1,3);
params.patch_transparency=0;
params.t_log_mode=0;
params.force_scale_factor=1;
params.fl_scale_factor=1;
params.force_record_ktr_offset=1;

% Update
params=parse_pv_pairs(params,varargin);

% Code

% Set output data
output.max_force=[];
output.min_force=[];
output.max_fl=[];
output.min_fl=[];
output.min_drawn_points=[];
output_drawn_indices=[];

% No of records
if (iscell(params.record_file_strings))
    no_of_records=numel(params.record_file_strings);
else
    no_of_records=1;
end

% Some updating
if (numel(params.t_offset_s)<no_of_records)
    params.t_offset_s=params.t_offset_s*ones(no_of_records,1);
end

for file_counter=1:no_of_records
    
    % Work out whether we are working with an slr file or
    % a slc file
    if (iscell(params.record_file_strings))
        file_string=params.record_file_strings{file_counter};
    else
        file_string=params.record_file_strings;
    end
    
    [~,~,file_extension]=fileparts(file_string);
    t=[];
    f=[];
    
    if (strcmp(file_extension,'.slr'))
        % It's an slr file
        d=ezread2(file_string);
        
        % Pull out the time and correct for time_offset
        if (numel(params.t_offset_s)>1)
            t=d.Time - params.t_offset_s(file_counter);
        else
            t=d.Time - params.t_offset_s;
        end
        
        % Pull off the force records by pCa
        for pCa_counter=1:length(params.slr_pCa_values)
            f_identifier = sprintf('pCa_%.0f_mean_force', ...
                10*params.slr_pCa_values(pCa_counter));
            f(pCa_counter,1:length(t))= ...
                params.force_scale_factor * d.(f_identifier);
            
            % Pull out the error range if necessary
            if (~isempty(params.force_error_label))
                f_error_identifier = sprintf( ...
                    'pCa_%.0f_%s_force', ...
                    10*params.slr_pCa_values(pCa_counter), ...
                    params.force_error_label);
                
                f_error(pCa_counter,1:length(t))= ...
                    params.force_scale_factor * d.(f_error_identifier);
            end
        end
        
        % Pull out the first fl_record
        fl_identifier=sprintf('pCa_%.0f_mean_fl', ...
            10*params.slr_pCa_values(pCa_counter));
        fl=d.(fl_identifier);
        
    else
        % It's an SLControl record
        d=transform_slcontrol_record( ...
            load_slcontrol_file(file_string), ...
            params.force_record_ktr_offset);
       
        % Pull out the time, force and the fl
        t=d.time;
        f=params.force_scale_factor * d.force';
        fl=params.fl_scale_factor * d.fl';
        
        pCa_value=d.pCa;
        
         % Correct for time_offset
        if (numel(params.t_offset_s)>1)
            t = t - params.t_offset_s(file_counter);
        else
            t = t - params.t_offset_s;
        end
        
        % There's no f_error
        f_error=[];
    end
        
    % Standard processing for both file types
    
    % Some defaults
    if (isempty(params.start_time_s))
        params.start_time_s=0;
    end
    if (isempty(params.stop_time_s))
        params.stop_time_s=t(end);
    end
    
    if (~params.t_log_mode)
        % This is the normal situation
        start_index=find(t>=params.start_time_s,1,'first');
        stop_index=find(t<=params.stop_time_s,1,'last');

        t_skip_points=1;
        t_indices=start_index:t_skip_points:stop_index;
        no_of_drawn_points=length(t_indices);

        while (no_of_drawn_points>=params.max_points_per_line)
            t_skip_points=t_skip_points+1;
            t_indices=start_index:t_skip_points:stop_index;
            no_of_drawn_points=length(t_indices);
        end
    else
        % This is the log mode
        start_index=find(t>(2*params.start_time_s),1,'first')
        stop_index=find(t<=params.stop_time_s,1,'last')
        
        t_skip_points=1;
        t_indices=start_index:t_skip_points:stop_index;
        no_of_drawn_points=length(t_indices);
        
        % Zero times to start
        t=t-params.start_time_s;
        
        % Log10 the t values you are going to plot
        t(t_indices)=log10(t(t_indices));
    end

    if (params.normalize_fl)
        fl=fl./fl(1);
    end
    
    % Draw the force record
    subplot(params.force_subplot);
    
    % Draw the error if required
    if (~isempty(params.force_error_label))
        for pCa_counter=1:length(params.slr_pCa_values)
            % Construct the f_box
            f_box = [f(pCa_counter,t_indices)+ ...
                            f_error(pCa_counter,t_indices) ...
                fliplr(f(pCa_counter,t_indices)- ...
                            f_error(pCa_counter,t_indices))]';
            t_box=[t(t_indices) ; flipud(t(t_indices))];
            
            % Draw the box
            patch(t_box,f_box, ...
                params.patch_color(file_counter,:), ...
                'EdgeColor','none');
            
            % Note the maximum
            output.max_force=max([output.max_force max(f_box)]);
            output.min_force=max([output.min_force min(f_box)]);
        end
    end
    
    % Plot
    [no_of_colors,~]=size(params.trace_colors);
    color_index=mod(file_counter,no_of_colors);
    if (color_index==0)
        color_index=no_of_colors;
    end
    plot(t(t_indices),f(:,t_indices),'-', ...
        'LineWidth',params.trace_line_width, ...
        'color',params.trace_colors(color_index,:));
    
    % Add in pCa values if appropriate
    if (params.display_pCa_values)
        if (strcmp(file_extension,'.slr'))
            % It's an slr file
            y_holder=[];
            for pCa_counter=1:length(params.slr_pCa_values)
                xt=t(t_indices(1))+params.pCa_x_offset* ...
                    (t(t_indices(end))-t(t_indices(1)));
                yt=f(pCa_counter,t_indices(end));
                y_holder=[y_holder yt];
                text(xt,yt, ...
                    sprintf('%.1f',params.slr_pCa_values(pCa_counter)), ...
                    'FontSize',params.pCa_font_size, ...
                    'HorizontalAlignment','center');
            end
            text(xt,params.pCa_label_y_offset*max(y_holder), ...
                'pCa','FontSize',params.pCa_font_size, ...
                'HorizontalAlignment','center');
        else
            % Check if there is a pCa value
            if (~isnan(pCa_value))
                % It's an slc file
                xt=t(t_indices(1))+params.pCa_x_offset* ...
                        (t(t_indices(end))-t(t_indices(1)));
                yt=f(1,t_indices(end));
                text(xt,yt, ...
                    sprintf('pCa %.1f',pCa_value), ...
                    'FontSize',params.pCa_font_size, ...
                    'HorizontalAlignment','center');
            end
        end
    end

    % Draw the fl record if required
    if (~isempty(params.fl_subplot))
        subplot(params.fl_subplot);
        
        if (isempty(params.single_fl_color))
           [no_of_colors,~]=size(params.trace_colors);
            color_index=mod(file_counter,no_of_colors);
            if (color_index==0)
                color_index=no_of_colors;
            end
            fl_color = params.trace_colors(color_index,:);
        else
            fl_color = params.single_fl_color;
        end

        plot(t(t_indices),fl(t_indices),'-', ...
            'LineWidth',params.trace_line_width, ...
            'color',fl_color);
    end
    
    % Output data
    output.last_force(file_counter) = f(:,t_indices(end));
    output.max_force=max([output.max_force max(f(:,t_indices))]);
    output.min_force=min([output.min_force min(f(:,t_indices))]);
    output.max_fl=max([output.max_fl max(fl(t_indices))]);
    output.min_fl=min([output.min_fl min(fl(t_indices))]);
    output.min_drawn_points=min([output.min_drawn_points ...
        min(t_indices)]);
end

% Limit axes
subplot(params.force_subplot);
if (~params.t_log_mode)
    xlim([t(start_index) t(stop_index)]);
else
%     xlim(log10([t(start_index) t(stop_index)]));
end
    

if (~isempty(params.fl_subplot))
    subplot(params.fl_subplot);
    xlim([t(start_index) t(stop_index)]);
end

% Update output data
output.min_time_s = t(start_index);
output.max_time_s = t(stop_index);
output.drawn_indices = t_indices;
        