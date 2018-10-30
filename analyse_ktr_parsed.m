function ktr_results=analyse_ktr_parsed(varargin)
% ktr_results=analyse_ktr(data_record,figure_number,max_fitting_period,
%                           log_mode,threshold,fit_double_exponential)
% Given an SLControl data record, display the ktr recovery in the figure
% specified by figure_number (default 1) with an x-axis scaled either linearly
% (log_mode=0, default) or log-wise (log_mode=1). Find the residual force
% (the minimum after the ktr re-stretch) and express it as a proportion of
% the steady-state isometric force. Find the maximum force in the
% max_fitting_period interval after re-stretch and thus deduce the relative
% ktr overshoot. Also find recovery rates, k_half (from the time to
% half_max force), k_tr (from a single exponential fit) and k_tr2 (from a
% double exponential fit). Double exponential fit is relatively time-consuming
% and is only performed if fit_double_exponential=1 (default 0). The relative residual
% has be less than recovery_threshold for a valid ktr calculation. The
% search for the residual is from the ktr peak (within 20 ms of re-stretch)
% for residual_search_factor * ktr_duration.

% Default values

% Default values

params.data=[];
params.figure_number=1;                 % figure number
params.max_fitting_period=10;           % max time span for curve fitting
params.log_scale=0;                     % 0 for linear time axis,
                                        % 1 for log
params.threshold=0.95;                  % deafult crossing time threshold
params.fit_double_exponential=0;        % 0 means don't fit (faster)
                                        % 1 means fit double exponential
params.residual_search_factor=1.5;      % number of ktr_step durations beyond
                                        % re-stretch to search for residual
params.recovery_threshold=1.0;          % residual has to be less than this
                                        % to be a valid ktr
params.spike_window=0.02;               % time window used to look for
                                        % spike height
params.draw_figure=1;                   % 0 means do not draw figure (quicker)
                                        % 1 means draw figure
params.steady_state_range=1:50;         % data points defining initial
                                        % steady_state
params.residual_search_smoothing_factor=0.001; 
                                       % multiple of sampling rate used to
                                       % smooth force record for residual
                                       % determination
params.max_smoothing_half_window=15;   % smoothing window used to establish max force


% Check for overrides
params=parse_pv_pairs(params,varargin);

% Now update function variables

data=params.data;
figure_number=params.figure_number;
max_fitting_period=params.max_fitting_period;
log_scale=params.log_scale;
threshold=params.threshold;
fit_double_exponential=params.fit_double_exponential;
residual_search_factor=params.residual_search_factor;
recovery_threshold=params.recovery_threshold;
spike_window=params.spike_window;
draw_figure=params.draw_figure;
steady_state_range=params.steady_state_range;
residual_search_smoothing_factor=params.residual_search_smoothing_factor;
max_smoothing_half_window=params.max_smoothing_half_window;

% Modification for ramps
if (data.ramp_mode==1)&(data.triangle_halftime<0.001)
    data.ktr_initiation_time=data.pre_triangle_time;
    data.ktr_duration=data.pre_ktr;
end
 
% Switch to appropriate figure

if (figure_number==0)
    draw_figure=0;
end
if (draw_figure)
    figure(figure_number);
    clf;
    hold on;
end

% Select data region to display

x_display_index_start=floor(data.sampling_rate*(data.ktr_initiation_time-data.ktr_duration));
x_display_index_stop=min([data.no_of_points ...
        ceil(data.sampling_rate*(data.ktr_initiation_time+2*data.ktr_duration+max_fitting_period))]);

if (log_scale==0)
    x_display_data=data.time(x_display_index_start:x_display_index_stop);
else
    x_display_data=log10([x_display_index_start:x_display_index_stop]-x_display_index_start+1);
end

if (draw_figure)
    plot(x_display_data,data.force(x_display_index_start:x_display_index_stop),'b');

    % Minimal text display
    [x_values]=get(gca,'XLim');
    [y_values]=get(gca,'YLim');
    text(x_values(1),y_values(1),truncate_display_string(data.file_string,50),'FontSize',8, ...
        'HorizontalAlignment','left','VerticalAlignment','bottom','Interpreter','none');
end

% Now find and display steady_state level

steady_state_level=mean(data.force(steady_state_range));

if (draw_figure)
    line([x_display_data(1) x_display_data(length(x_display_data))],[steady_state_level steady_state_level], ...
        'Color',[1 0 1]);
    text(x_display_data(end),steady_state_level,'Steady-state force', ...
        'Color',[1 0 1]);
end

% Now find and display base-line level

baseline_range=[ceil(data.sampling_rate*(data.ktr_initiation_time+0.5*data.ktr_duration)): ...
        floor(data.sampling_rate*(data.ktr_initiation_time+data.ktr_duration))-1];
baseline_level=mean(data.force(baseline_range));

if (draw_figure)
    line([x_display_data(1) x_display_data(length(x_display_data))],[baseline_level baseline_level], ...
        'Color',[1 1 0]);
    text(x_display_data(end),baseline_level,'Baseline force', ...
        'Color',[1 1 0]);
end

% Now find residual point by first finding the maximum after the ktr re-stretch
% and then looking for the subsequent minimum

resid_search_min_index=floor(data.sampling_rate*(data.ktr_initiation_time+0.5*data.ktr_duration));
resid_search_max_index=ceil(data.sampling_rate* ...
    (data.ktr_initiation_time+(0.5+residual_search_factor)*data.ktr_duration));
resid_search_range=[resid_search_min_index:resid_search_max_index];
resid_extra=ceil(residual_search_smoothing_factor*data.sampling_rate);
resid_search_force= ...
    smooth_array_in_place(data.force(resid_search_min_index-resid_extra:resid_search_max_index+resid_extra), ...
    resid_extra);
resid_search_force=resid_search_force(resid_extra+1:resid_extra+length(resid_search_range));

% Display search_range

if (draw_figure)
    if (log_scale==0)
        line(data.time(resid_search_range),resid_search_force,'Color',[1 0 0]);
    else
        line(log10((resid_search_range)-x_display_index_start+1),resid_search_force,'Color',[1 0 0]);
    end
end

% Find residual

% Start by searching for the last minimum
pks = peakfinder( ...
    resid_search_force, ...
    0.01*(max(resid_search_force)-min(resid_search_force)), ...
    max(resid_search_force),-1,true,false);

% Check if you didn't find one.
if ((isempty(pks))|| ...
        ((numel(pks)==1)&&(resid_search_force(pks(1))==min(resid_search_force))))
    dy = diff(resid_search_force,1);
    [~,resid_index]=max(dy);
else
    resid_index = pks(end);
end
resid_index = resid_search_min_index + resid_index - 1;
resid_time = data.time(resid_index);
residual_y_value = data.force(resid_index);
restretch_spike_height = max(resid_search_force);

% Display residual

if (draw_figure)
    line([x_display_data(1) x_display_data(length(x_display_data))],residual_y_value*ones(1,2),'Color',[ 0 1 0]);
    if (log_scale==0)
        line(data.time(resid_index)*ones(1,2),[baseline_level residual_y_value],'Color',[0 1 0]);
    else
        line(log10(resid_index-x_display_index_start+1)*ones(1,2), ...
            [baseline_level residual_y_value],'Color',[0 1 0]);
        plot(log10(resid_index-x_display_index_start+1),residual_y_value,'go');
    end
end

% Find max_level

max_search_range=[resid_index:x_display_index_stop];
[max_y_value,index]=max(smooth_array_in_place(data.force(max_search_range),max_smoothing_half_window));
max_index=resid_index+index-1;
max_time=data.time(max_index);

% Display max_level

if (draw_figure)
    line([x_display_data(1) x_display_data(length(x_display_data))],max_y_value*ones(1,2),'Color',[1 1 0]);
    if (log_scale==0)
        line(data.time(max_index)*ones(1,2),[baseline_level max_y_value],'Color',[0 1 0]);
    else
        line(log10(max_index-x_display_index_start+1)*ones(1,2), ...
            [baseline_level max_y_value],'Color',[0 1 0]);
    end
end

% Update ktr_results

ktr_results.resid_time=(resid_time-(data.ktr_initiation_time+data.ktr_duration));
ktr_results.max_time=(max_time-(data.ktr_initiation_time+data.ktr_duration));
ktr_results.recovery_time=max_time-resid_time;

% Extract ktr_data for fitting

x_fit_data=data.time(resid_index:max_index)-data.time(resid_index);
y_fit_data=data.force(resid_index:max_index)-data.force(resid_index);
y_fit_scaling_factor=max([max(y_fit_data) 1]);        % prevents problem with decaying trace
y_fit_data=y_fit_data./y_fit_scaling_factor;

% Calculate and draw k_half

half_point=max([min(find(y_fit_data>0.5)) 2]);        % prevents problem with decaying trace  
half_point_index=half_point+resid_index-1;
half_time=data.time(half_point_index)-data.time(resid_index);
k_half=-log(0.5)/half_time;
half_point_y_value=0.5*(max_y_value+residual_y_value);

if (draw_figure)
    line([x_display_data(1) x_display_data(length(x_display_data))],half_point_y_value*ones(1,2), ...
        'Color',[1 1 0]);
    if (log_scale==0)
        line(data.time(half_point_index)*ones(1,2),[baseline_level half_point_y_value],'Color',[0 1 0]);
    else
        line(log10(half_point_index-x_display_index_start+1)*ones(1,2), ...
            [baseline_level half_point_y_value],'Color',[0 1 0]);
    end
end

% Calculate single exponential fit

x_ktr=x_fit_data;
y_ktr=data.force(resid_index:max_index)-data.force(resid_index);

p(1)=y_fit_scaling_factor;
p(2)=k_half;

p=fminsearch(@single_exponential_recovery_fit,p,[],x_ktr,y_ktr);
y_exponential=p(1)*(1-exp(-p(2).*x_ktr));

if (draw_figure)
    if (log_scale==0)
        line(data.time(resid_index:max_index),residual_y_value+y_exponential,'Color',[1 0 1]);
    else
        line(log10([resid_index:max_index]-x_display_index_start+1)',residual_y_value+y_exponential,'Color',[1 0 1]);
    end
end

% Update ktr_results

ktr_results.single_exp_data.p=p;
ktr_results.single_exp_data.r_squared=calculate_r_squared(y_exponential,y_ktr);

% Calculate double exponential fit

if (fit_double_exponential)
	x_ktr2=x_fit_data;
	y_ktr2=data.force(resid_index:max_index)-data.force(resid_index);
	
	p(1)=0.75*y_fit_scaling_factor;
	p(2)=0.5*k_half;
	p(3)=0.25*y_fit_scaling_factor;
	p(4)=2*k_half;
	
	p=fminsearch(@double_exponential_recovery_fit,p,[],x_ktr2,y_ktr2);
	y_doub_exponential=p(1)*(1-exp(-p(2).*x_ktr2))+p(3)*(1-exp(-p(4).*x_ktr2));

    if (draw_figure)
        if (log_scale==0)
            line(data.time(resid_index:max_index),residual_y_value+y_doub_exponential,'Color',[0 1 1]);
        else
            line(log10([resid_index:max_index]-x_display_index_start+1)',residual_y_value+y_doub_exponential,'Color',[0 1 1]);
        end
    end
    
	% Update ktr_results
	
	ktr_results.double_exp_data.p=p;
	ktr_results.double_exp_data.r_squared=calculate_r_squared(y_doub_exponential,y_ktr2);
end

% Calculate and draw crossing_point

crossing_point=min(find((data.force(resid_index:max_index)-data.force(resid_index))> ...
    (threshold*(steady_state_level-data.force(resid_index)))));
crossing_point_index=crossing_point+resid_index-1;
crossing_point_time=data.time(crossing_point_index)-data.time(resid_index);
crossing_point_y_value=data.force(crossing_point_index);

if (draw_figure)
    if (~isnan(crossing_point_y_value))
        line([x_display_data(1) x_display_data(length(x_display_data))],crossing_point_y_value*ones(1,2), ...
            'Color',[1 1 0]);
        if (log_scale==0)
            line(data.time(crossing_point_index)*ones(1,2),[baseline_level crossing_point_y_value],'Color',[0 1 0]);
            text(data.time(crossing_point_index),crossing_point_y_value,sprintf('%.3f',threshold), ...
                'HorizontalAlignment','left','VerticalAlignment','top','FontSize',6);
        else
            line(log10(crossing_point_index-x_display_index_start+1)*ones(1,2), ...
                [baseline_level crossing_point_y_value],'Color',[0 1 0]);
            text(log10(crossing_point_index-x_display_index_start+1),crossing_point_y_value,sprintf('%.3f',threshold), ...
                'HorizontalAlignment','left','VerticalAlignment','top','FontSize',6);
        end
    end
end

% Update ktr_results

% Filter ktr_results to NaN if resid > recovery_threshold

ktr_results.isometric_tension=steady_state_level-baseline_level;
ktr_results.residual=(residual_y_value-baseline_level)/ktr_results.isometric_tension;
ktr_results.relative_restretch_spike_height=(restretch_spike_height-baseline_level)/ktr_results.isometric_tension;

if (length(recovery_threshold)==1)
    if (ktr_results.residual<=recovery_threshold)& ...
            ((k_half<(data.sampling_rate/20))&(ktr_results.single_exp_data.p(2)<(data.sampling_rate/20)))
        ktr_results.k_half=k_half;
        ktr_results.k_tr=ktr_results.single_exp_data.p(2);
        ktr_results.crossing_point_time=crossing_point_time;
        ktr_results.overshoot=(max_y_value-baseline_level)/ktr_results.isometric_tension;
    else
        ktr_results.k_half=NaN;
        ktr_results.k_tr=NaN;
        ktr_results.crossing_point_time=NaN;
        ktr_results.overshoot=NaN;
    end
else
    if ((ktr_results.residual>recovery_threshold(1))&(ktr_results.residual<=recovery_threshold(2))) ...
            &((k_half<(data.sampling_rate/20))&(ktr_results.single_exp_data.p(2)<(data.sampling_rate/20)))
        ktr_results.k_half=k_half;
        ktr_results.k_tr=ktr_results.single_exp_data.p(2);
        ktr_results.crossing_point_time=crossing_point_time;
        ktr_results.overshoot=(max_y_value-baseline_level)/ktr_results.isometric_tension;
    else
        ktr_results.k_half=NaN;
        ktr_results.k_tr=NaN;
        ktr_results.crossing_point_time=NaN;
        ktr_results.overshoot=NaN;
    end
end 


% Update display

if (draw_figure)
    drawnow;
end
