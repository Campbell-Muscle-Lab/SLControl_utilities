function transformed_record=transform_slcontrol_record(varargin)
% transformed_record=transform_slcontrol_record(data)
% Duplicates SLControl->AnalysisDisplay->Transform->Change to
% Calibrated Values
% Default operation changes force to N m^-2 and zeros force relative
% to the points 60 to 90% of the way through the ktr step.
% Expresses SL in metres (assuming the first point is the resting
% value) and changes FL to metres.
% If the second argument is 'no_zero' the function changes force
% to N m^-2 but doesn't zero the record

no_of_arguments=length(varargin);

switch no_of_arguments
    case {1}
        data=varargin{1}(:);
        offset_force=1;
    case {2}
        data=varargin{1}(:);
        offset_force=varargin{2}(:);
    otherwise
        error('Invalid number (%i) of arguments for transform_slcontrol_record', ...
            no_of_arguments);
end

% Transform force record

if (offset_force==1)
    % Now deduce the ktr_range
    if (isnan(data.ktr_initiation_time))
        data.ktr_initiation_time = data.hold_time + data.pre_triangle_time + (2*data.no_of_triangles*data.triangle_halftime) + ...
            ((data.no_of_triangles)*data.inter_triangle_interval) + data.pre_ktr+ .01
    end
    
    %re-calculate sampling rate
    zero_sec= find(data.time == 0)
    one_sec= find(data.time == 1.0000)
    data.sampling_rate= (one_sec-zero_sec)
    
%     ktr_begin_point=floor(data.sampling_rate* ...
%         (data.ktr_initiation_time+(0.6*data.ktr_duration)))
%     ktr_end_point=ceil(data.sampling_rate* ...
%         (data.ktr_initiation_time+(0.9*data.ktr_duration)))
    velocity= diff(data.fl)
%     fc_plot= subplot(3,1,1)
%     plot(velocity)
    [~,ktr_start] =max(velocity)
    disp(ktr_start)
    ktr_begin_point=floor(data.sampling_rate*(0.6*data.ktr_duration) + ktr_start)
    ktr_end_point=ceil(data.sampling_rate*(0.9*data.ktr_duration) + ktr_start)
    
    mean_force=mean(data.force(ktr_begin_point:ktr_end_point))
    data.force=(data.force-mean(data.force(ktr_begin_point:ktr_end_point))).* ...
        (data.FORCE_CALIBRATION/data.area);
    
%     %Normalize Force to intitial Value
%     initial_force=data.force(1)
%     data.force_normalized= (data.force)/(initial_force);

    %plot force trace with ktr_points
%     subplot(3,1,2)
%     plot(data.force)
%     hold on
%     plot(ktr_begin_point,data.force(ktr_begin_point),'r*')
%     plot(ktr_end_point,data.force(ktr_end_point),'b*')
%     hold off
    %zoom to ktr points
%     ktr_begin_point_padding= ktr_begin_point - (data.sampling_rate* 0.075)
%     ktr_end_point_padding= ktr_end_point + (data.sampling_rate* 0.075)
%     subplot(3,1,3)
%     plot(data.time(ktr_begin_point_padding:ktr_end_point_padding),data.force(ktr_begin_point_padding:ktr_end_point_padding))
%     hold on
%     plot(data.time(ktr_begin_point),data.force(ktr_begin_point),'r*')
%     plot(data.time(ktr_end_point),data.force(ktr_end_point),'b*')
%     hold off
    
%     savefig
    %create output file for figure
%     split_file_name= split(data.file_string,'\')
%     file_location= strjoin(split_file_name(1:end-1),{'\'})
%     prep_name= string(split_file_name(end-1))
%     pCa_string_figure= strrep(string(data.pCa), '.','_')
%     filename_string= strcat(string(file_location), '\', prep_name, '_', pCa_string_figure, '_force_baseline_correction')
%     figure_export('output_file_string', string(filename_string),'output_type','png')
    
else
    if (offset_force == 0)
        data.force=data.force.*(data.FORCE_CALIBRATION/data.area);
    end
    if (offset_force == -1)
        data.force = data.force.*(data.FORCE_CALIBRATION/data.area);
        data.force = data.force-min(data.force);
    end
end
    
% Offset and scale SL if appropriate
if (isfield(data,'sl'))
    fractional_fl_change=(data.sl-data.sl(1)).*data.FL_COMMAND.*data.FL_POLARITY./ ...
        (data.muscle_length*data.SL_CALIBRATION);
    data.sl=(1+fractional_fl_change)*data.sarcomere_length;
end

% Offset and scale FL

data.fl=((data.fl-data.fl(1)).*data.FL_POLARITY.*data.FL_RESPONSE)+data.muscle_length;

% %Calculate length relative to initial length
for i= 1:numel(data.fl)
data.length_change(i,1)= data.fl(i)/data.fl(1);
end
% Now copy t to transformed_record

transformed_record=data
