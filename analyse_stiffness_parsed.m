function stiffness_results=analyse_stiffness_parsed(varargin);
% ktr_results=analyse_stiffness(data_record,figure_number)
% Given an SLControl data record, display the ktr recovery in the figure
% specified by figure_number (default 1)

% Default values
% Check and set input arguments

params.data=[];
params.figure_number=1;
params.breakpoint_figure=2;
params.padding_time=0.075;
params.initial_fitting_time=0.075;
params.second_fitting_time=0.1;
params.SREC_detection_threshold=2;

params.smoothing_window=15;
params.base_smoothing_points=25;

params.apply_smoothing_for_regression = 1;

% Update data
params=parse_pv_pairs(params,varargin)

% Set function variables
data=params.data;
figure_number=params.figure_number;
padding_time=params.padding_time;
initial_fitting_time=params.initial_fitting_time;
second_fitting_time=params.second_fitting_time;
SREC_detection_threshold=params.SREC_detection_threshold;

smoothing_window=params.smoothing_window;
base_smoothing_points=params.base_smoothing_points;

% Check whether SL_CALIBRATION is valid
if (abs(data.SL_CALIBRATION-999.9)>eps)
    valid_sl_calibration=1;
else
    valid_sl_calibration=0;
end

% Options

optimise_options=optimset('MaxIter',200,'MaxFunEvals',2500);

% Switch to appropriate figure

figure(figure_number);
clf;

% Set triangle start, stop and peak times and indices

tri_start_time(1)=data.pre_triangle_time;
for tri_counter=2:data.no_of_triangles
    if (tri_counter==2)
        tri_start_time(tri_counter)=tri_start_time(1)+ ...
            ((2-data.ramp_mode)*data.relative_first_triangle_size*data.triangle_halftime)+ ...
            data.inter_triangle_interval;
    else
        tri_start_time(tri_counter)=tri_start_time(tri_counter-1)+ ...
            ((2-data.ramp_mode)*data.triangle_halftime)+data.inter_triangle_interval;
    end
end

tri_peak_time(1)=tri_start_time(1)+data.relative_first_triangle_size*data.triangle_halftime;
for tri_counter=2:data.no_of_triangles
    tri_peak_time(tri_counter)=tri_start_time(tri_counter)+data.triangle_halftime;
end

tri_stop_time(1)=tri_start_time(1)+ ...
    (2-data.ramp_mode)*data.relative_first_triangle_size*data.triangle_halftime;
for tri_counter=2:data.no_of_triangles
    tri_stop_time(tri_counter)=tri_start_time(tri_counter)+(2-data.ramp_mode)*data.triangle_halftime;
end

x_display_index_start=max([1 floor(data.sampling_rate*(tri_start_time(1)-padding_time))]);
x_display_index_stop=min([floor(data.sampling_rate*(tri_stop_time(data.no_of_triangles)+padding_time))]);


% Display with minimal text
subplot(3,2,1);
plot(data.time(x_display_index_start:x_display_index_stop),data.force(x_display_index_start:x_display_index_stop),'b-');
xlabel('Time (s)');
ylabel('Force');

subplot(3,2,2);
plot(data.sl(x_display_index_start:x_display_index_stop),data.force(x_display_index_start:x_display_index_stop),'b-');
xlabel('Force');
ylabel('SL');

subplot(3,2,3);
plot(data.time(x_display_index_start:x_display_index_stop),data.sl(x_display_index_start:x_display_index_stop),'b-');
xlabel('Time (s)');
ylabel('SL');

subplot(3,2,4);
plot(data.fl(x_display_index_start:x_display_index_stop),data.force(x_display_index_start:x_display_index_stop),'b-');
xlabel('FL');
ylabel('Force');

subplot(3,2,5);
plot(data.time(x_display_index_start:x_display_index_stop),data.fl(x_display_index_start:x_display_index_stop),'b-');
xlabel('Time (s)');
ylabel('FL');
[x_values]=get(gca,'XLim');
[y_values]=get(gca,'YLim');
text(x_values(1),y_values(1),truncate_display_string(data.file_string,35),'FontSize',8, ...
    'HorizontalAlignment','left','VerticalAlignment','bottom','Interpreter','none');

subplot(3,2,6);
plot(data.fl(x_display_index_start:x_display_index_stop),data.sl(x_display_index_start:x_display_index_stop),'b-');
xlabel('SL');
ylabel('FL');

% Analysis

for tri_counter=1:data.no_of_triangles
    
    % Peak values
    peak_index=ceil(data.sampling_rate*tri_peak_time(tri_counter));
    stiffness_results(tri_counter,1,1).peak_tri_tension=data.force(peak_index);
    subplot(3,2,1);
    line([tri_peak_time(tri_counter)-0.05 tri_peak_time(tri_counter)+0.05], ...
        stiffness_results(tri_counter,1,1).peak_tri_tension*ones(2,1),'Color',[0 1 0]);
    
    % P_end - force value at end of triagnle
    end_index=ceil(data.sampling_rate*tri_stop_time(tri_counter));
    stiffness_results(tri_counter,1,1).P_end = data.force(end_index);
    
    for phase_counter=1:2
        % Stiffness calculations
        for mode_counter=1:2
            
            % Set indices
            if (phase_counter==1)
                % Initial stiffness
                indices=[floor(data.sampling_rate*tri_start_time(tri_counter)): ...
                ceil(data.sampling_rate*(tri_start_time(tri_counter)+initial_fitting_time))];
            else
                % Second stiffness
                indices=[floor(data.sampling_rate*(tri_peak_time(tri_counter)-second_fitting_time)): ...
                ceil(data.sampling_rate*(tri_peak_time(tri_counter)))];
            end
            
            % FL or SL calculation
            if (mode_counter==1)
                % Force/FL regression
                x_data=data.fl(indices);
                y_data=data.force(indices);
                if (params.apply_smoothing_for_regression)
                    x_data = smooth_array_in_place(x_data,smoothing_window);
                    y_data = smooth_array_in_place(y_data,smoothing_window);
                end
                prep_length=data.muscle_length;
                subplot(3,2,4);
            else
                % Force/SL regression
                if (valid_sl_calibration)
                    x_data=data.sl(indices);
                    y_data=data.force(indices);
                    prep_length=data.sarcomere_length;
                    subplot(3,2,2);
                end
            end
            
            if ((mode_counter==1) | ...
                    ((mode_counter==2)&(valid_sl_calibration)))

                % Run Regression analysis
                [gradient,intercept,r] = fit_robust_line(x_data,y_data);

                % Draw fit
                line(x_data,y_data,'Color',[0 1 0]);
                line([x_data(1) x_data(end)], ...
                    intercept+gradient*[x_data(1) x_data(end)], ...
                    'Color',[1 0 0]);
            
                % Update results
                
                if (phase_counter==1)
                    if (mode_counter==1)
                        % FL
                        stiffness_results(tri_counter).SREC_YM_fl= ...
                            gradient*data.muscle_length;
                        stiffness_results(tri_counter).SREC_stiff_fl= ...
                            gradient;
                    else
                        % SL
                        stiffness_results(tri_counter).SREC_YM_sl= ...
                            gradient*data.sarcomere_length;
                        stiffness_results(tri_counter).SREC_stiff_sl= ...
                            gradient;
                    end
                else
                    % Final Phase
                    if (mode_counter==1)
                        % FL
                        stiffness_results(tri_counter).final_YM_fl = ...
                            gradient*data.muscle_length;
                        stiffness_results(tri_counter).final_stiff_fl = ...
                            gradient;
                    else
                        % SL
                        stiffness_results(tri_counter).final_YM_sl = ...
                            gradient*data.sarcomere_length;
                        stiffness_results(tri_counter).final_stiff_sl = ...
                            gradient;
                    end
                end
            else
                % No SL calibration
                stiffness_results(tri_counter).SREC_YM_sl = NaN;
                stiffness_results(tri_counter).SREC_stiff_sl = NaN;
                stiffness_results(tri_counter).final_YM_sl = NaN;
                stiffness_results(tri_counter).final_stiff_sl = NaN;
            end
        end
    end
end

% SREC analysis

% Switch to appropriate figure
breakpoint_figure=figure(params.breakpoint_figure);
clf;

% Set up display
no_of_rows=3;
no_of_columns=2*data.no_of_triangles;

for tri_counter=1:data.no_of_triangles
   
    % Set display indices
    x_display_index_start=max([1 floor(data.sampling_rate* ...
        (tri_start_time(tri_counter)-padding_time))]);
    x_display_index_stop=min([floor(data.sampling_rate* ...
        (tri_peak_time(tri_counter)+padding_time)) data.no_of_points]);

    % Calculate t,sm_force and sm_dfdt
    t=data.time(x_display_index_start:x_display_index_stop);
    f=data.force(x_display_index_start:x_display_index_stop);
    fl=data.fl(x_display_index_start:x_display_index_stop);
    sl=data.sl(x_display_index_start:x_display_index_stop);
    
    t_range=t(end)-t(1);
    t_padding=0.1*t_range;

    sm_f=(smooth_array_in_place(data.force(x_display_index_start:x_display_index_stop),smoothing_window))';

    % Display
    subplot(no_of_rows,no_of_columns,2*tri_counter-1);
    plot(t,f,'b',t,sm_f,'m');
    xlim([t(1)-t_padding t(end)+t_padding]);
    hold on;

    subplot(no_of_rows,no_of_columns,2*tri_counter);
    hold on;

    subplot(no_of_rows,no_of_columns,no_of_columns+2*tri_counter-1);
    plot(t,sl,'b');
    xlim([t(1)-t_padding t(end)+t_padding]);
    hold on;

    subplot(no_of_rows,no_of_columns,no_of_columns+2*tri_counter);
    plot(t,fl,'b');
    xlim([t(1)-t_padding t(end)+t_padding]);
    hold on;

    % Start analysis

    SREC_point=0;

    % Set indices
    tri_start_index=ceil(data.sampling_rate*padding_time);
    if (tri_counter==1)
        tri_stop_index=ceil(data.sampling_rate*(padding_time+ ...
            data.relative_first_triangle_size*data.triangle_halftime));
    else
        tri_stop_index=ceil(data.sampling_rate* ...
            (padding_time+data.triangle_halftime));
    end

    % Display
    subplot(no_of_rows,no_of_columns,no_of_columns+2*tri_counter);
    plot([t(tri_start_index) t(tri_stop_index)],[fl(tri_start_index) fl(tri_stop_index)],'m^');
    
    % Look for a peak by subtracting an estimate of the passive force
    % and seeing if the remaining cb_force exceeds the noise threshold
    
    % First take off final slope
    base_fl=mean(fl(tri_start_index - ...
        base_smoothing_points:tri_start_index));
    base_force=mean(f(tri_start_index - ...
        base_smoothing_points:tri_start_index));

    sm_fl=smooth_array_in_place(fl,smoothing_window);
    sm_fl(1:smoothing_window)=sm_fl(smoothing_window+1);
    passive_force=base_force + ...
        (sm_fl-base_fl)*stiffness_results(tri_counter).final_stiff_fl;
    
    % Now calculate cb_stretch_force
    cb_stretch_force=f-passive_force;

    subplot(no_of_rows,no_of_columns,2*tri_counter-1);
    plot(t,passive_force,'g')

    subplot(no_of_rows,no_of_columns,2*tri_counter);
    plot(t,cb_stretch_force,'b');
    sm_cb_stretch_force=smooth_array_in_place(cb_stretch_force, ...
        smoothing_window);
    plot(t,sm_cb_stretch_force,'m');
 
    noise_level=mean(abs(cb_stretch_force-sm_cb_stretch_force));
    plot([t(1) t(end)],[noise_level noise_level],'g');
    plot([t(1) t(end)], ...
        SREC_detection_threshold*[noise_level noise_level],'c');

    [max_cb_stretch_force,max_cb_stretch_force_point]= ...
        max(sm_cb_stretch_force);

    if (max_cb_stretch_force>SREC_detection_threshold*noise_level)
        
        % Valid SREC detected

        SREC_point=max_cb_stretch_force_point;

        base_force=mean(sm_f(tri_start_index-base_smoothing_points:tri_start_index));
        base_fl=mean(fl(tri_start_index-base_smoothing_points:tri_start_index));
        base_sl=mean(sl(tri_start_index-base_smoothing_points:tri_start_index));
        
        % Update graphics
        
        subplot(no_of_rows,no_of_columns,2*tri_counter-1);
        text(t(tri_start_index),sm_f(tri_stop_index), ...
            'SREC detected');
        plot([t(SREC_point)],[sm_f(SREC_point)],'mv');
        
        subplot(no_of_rows,no_of_columns,2*tri_counter);
        plot([t(SREC_point)],[sm_cb_stretch_force(SREC_point)],'mv');
        
        if (valid_sl_calibration)
            subplot(no_of_rows,no_of_columns,no_of_columns+2*tri_counter-1);
            plot([t(SREC_point)],[sl(SREC_point)],'mv');
        end

        subplot(no_of_rows,no_of_columns,no_of_columns+2*tri_counter);
        plot([t(SREC_point)],[fl(SREC_point)],'mv');
        
        % Update stiffness_results
        
        stiffness_results(tri_counter).SREC_f = ...
            sm_f(SREC_point)-base_force;
        stiffness_results(tri_counter).SREC_fl = ...
            stiffness_results(tri_counter).SREC_f/ ...
                stiffness_results(tri_counter).SREC_stiff_fl;
        if (valid_sl_calibration)
            stiffness_results(tri_counter).SREC_sl = ...
                stiffness_results(tri_counter).SREC_f/ ...
                    stiffness_results(tri_counter).SREC_stiff_sl;
        else
            stiffness_results(tri_counter,1,1).SREC_sl=NaN;
        end

        % Update graphics
        
        if (valid_sl_calibration)
            subplot(no_of_rows,no_of_columns,2*no_of_columns+2*tri_counter-1);
            hold on;
            plot(sl,f,'b',sl,sm_f,'m')
            plot([min(sl) max(sl)],[base_force base_force],'g-');
            plot([min(sl) max(sl)],[sm_f(SREC_point) sm_f(SREC_point)],'g-');
            plot([min(sl) min(sl)+stiffness_results(tri_counter).SREC_sl], ...
                [base_force sm_f(SREC_point)],'g-');
        end

        subplot(no_of_rows,no_of_columns,2*no_of_columns+2*tri_counter);
        hold on;
        plot(fl,f,'b',fl,sm_f,'m')
        plot([min(fl) max(fl)],[base_force base_force],'g-');
        plot([min(fl) max(fl)],[sm_f(SREC_point) sm_f(SREC_point)],'g-');
        plot([min(fl) min(fl)+stiffness_results(tri_counter).SREC_fl], ...
            [base_force sm_f(SREC_point)],'g-');

    else
                
        % SREC not found
        
        subplot(no_of_rows,no_of_columns,2*tri_counter-1);
        text(t(tri_start_index),sm_f(tri_stop_index),'No peak detected');

        stiffness_results(tri_counter).SREC_f=NaN;
        stiffness_results(tri_counter).SREC_sl=NaN;
        stiffness_results(tri_counter).SREC_fl=NaN;

    end
end