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
    ktr_begin_point=floor(data.sampling_rate* ...
        (data.ktr_initiation_time+(0.6*data.ktr_duration)));
    ktr_end_point=ceil(data.sampling_rate* ...
        (data.ktr_initiation_time+(0.9*data.ktr_duration)));
    % Offset and scale the force
    data.force=(data.force-mean(data.force(ktr_begin_point:ktr_end_point))).* ...
        (data.FORCE_CALIBRATION/data.area);
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

% Now copy t to transformed_record

transformed_record=data;
