function display_slcontrol_record(varargin)
% display_sl_control_record(data,figure_number,log_mode,text_display,[channels])
% Given a data record (data), display_slcontrol_record shows an SLControl
% record in figure number (figure_number, default=1) with either linear
% (log_mode=0, default) or log (log_mode=1) x-axis scaling. text_display=1 (default) outputs
% some descriptive data about the file. Only channel numbers listed in [channels] (1-3) are
% displayed (default [1 2 3]). Only the data record (data) must be supplied. All other
% parameters are optional

% Text formatting variables

no_of_text_lines=8;
max_line_length=40;

% Check and set input arguments

no_of_arguments=length(varargin);

switch no_of_arguments
    case {1}
        data=varargin{1}(1);
        figure_number=1;
        log_mode=0;
        text_display=1;
        channels=[1 2 3];
    case {2}
        data=varargin{1}(1);
        figure_number=varargin{2}(1);
        log_mode=0;
        text_display=1;
        channels=[1 2 3];
    case {3}
        data=varargin{1}(1);
        figure_number=varargin{2}(1);
        log_mode=varargin{3}(1);
        text_display=1;
        channels=[1 2 3];
    case {4}
        data=varargin{1}(1);
        figure_number=varargin{2}(1);
        log_mode=varargin{3}(1);
        text_display=varargin{4}(1);
        channels=[1 2 3];
    case {5}
        data=varargin{1}(1);
        figure_number=varargin{2}(1);
        log_mode=varargin{3}(1);
        text_display=varargin{4}(1);
        channels=varargin{5};
    otherwise
        error('Invalid number (%i) of arguments for display_slcontrol_record()',no_of_arguments);
end
no_of_channels=length(channels);

% Stop and return if the figure_number==0
if (figure_number==0)
    return;
end

figure(figure_number);
clf;

% Display channels in turn
% Code only implemented for 1 to 3 channels at this time
% Set 

for counter=1:no_of_channels
    sp=subplot(no_of_channels,1,counter);
    
    if (log_mode==0)
        x_data=data.time;
    else
        x_data=log10(1:data.no_of_points);
    end
           
    y_data=return_y_data(data,channels(counter));
        
    line(x_data,y_data,'Color',[0 0 1]);
    ylabel(return_y_title(channels(counter)));
end

% Title x-axis

subplot(no_of_channels,1,no_of_channels);
if (log_mode==0)
    xlabel('Time (s)');
else
    xlabel('log10(Time(s))');
end

% Display diagnostic text if appropriate

if (text_display==1)&(isfield(data,'version'))
    
    subplot(no_of_channels,1,no_of_channels);
    [x_values]=get(gca,'XLim');
    [y_values]=get(gca,'YLim');
    
    y_top=y_values(2);
    y_spacing=(y_values(2)-y_values(1))/(no_of_text_lines-1);
    y_position=y_top;
    
    output_string=sprintf('%s',truncate_display_string(data.file_string,max_line_length));
    text(mean(x_values),y_position,output_string,'Interpreter','none','FontSize', 8, ...
        'HorizontalAlignment','left','VerticalAlignment','top');
    y_position=y_position-y_spacing;
    if (length(data.file_info_string)>0)
        text(mean(x_values),y_position, ...
            truncate_display_string(data.file_info_string,max_line_length), ...
            'Interpreter','none','FontSize', 8, ...
            'HorizontalAlignment','left','VerticalAlignment','top');
        y_position=y_position-y_spacing;                        
    end
    output_string=sprintf('%i:%i:%i   %i/%i/%i', ...
        data.creation_time.hours, ...
        data.creation_time.mins, ...
        data.creation_time.secs, ...
        data.creation_time.day, ...
        data.creation_time.month, ...
        data.creation_time.year);
    text(mean(x_values),y_position,output_string,'Interpreter','none','FontSize', 8, ...
        'HorizontalAlignment','left','VerticalAlignment','top');
    y_position=y_position-y_spacing;
    output_string=strrep(sprintf('FL (m): %4.3e   SL (m): %g', ...
        data.muscle_length,data.sarcomere_length),'e-00','e-');
    text(mean(x_values),y_position,output_string,'Interpreter','none','FontSize', 8, ...
        'HorizontalAlignment','left','VerticalAlignment','top');
    y_position=y_position-y_spacing;
    output_string=strrep(sprintf('Area (m^2): %g',data.area),'e-00','e-');
    text(mean(x_values),y_position,output_string,'Interpreter','none','FontSize', 8, ...
        'HorizontalAlignment','left','VerticalAlignment','top');
    y_position=y_position-y_spacing;
    output_string=strrep(sprintf('COM (um V^-1): %4.3e  RES: %4.3e', ...
        data.FL_COMMAND,data.FL_RESPONSE),'e-00','e-');
    text(mean(x_values),y_position,output_string,'Interpreter','none','FontSize', 8, ...
        'HorizontalAlignment','left','VerticalAlignment','top');

end

% Update display

drawnow;
    
% Sub-functions below

function y_data=return_y_data(data,channel_number);
% Code returns appropriate data channel for the channel number

switch channel_number;
    case {1}
        y_data=data.force;
    case {2}
        y_data=data.sl;
    case {3}
        y_data=data.fl;
    case {4}
        y_data=data.intensity;
    case {5}
        y_data=data.stimulus;
    otherwise
        y_data=data.time;
end

function y_title=return_y_title(channel_number);
% Code returns appropriate y axis title

switch channel_number;
    case {1}
        y_title=sprintf('Force');
    case {2}
        y_title=sprintf('Sarcomere Length');
    case {3}
        y_title=sprintf('Fibre Length');
    case {4}
        y_title=sprintf('Intensity');
    case {5}
        y_title=sprintf('Stimulus');
    otherwise
        y_title=sprintf('Default');
end


