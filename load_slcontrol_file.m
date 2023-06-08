function data_file=load_slcontrol_file(varargin)
% data_file=load_slcontrol_file(file_string,header_only)
% Extracts channel records and most of the header file information from the
% SLCONTROL file specified by file_string
% The file_string should include the path if the data file is not in the local directory)
% Data is returned as a structure - lengths are expressed as meters,
% times as seconds, creation time and channel data are returned as
% structures within the data_file structure itself
% If header_only==1, function only scans the header
% Please contact Ken Campbell (k.s.campbell@uky.edu) if you need further help or
% information

no_of_arguments=length(varargin);

switch no_of_arguments
    case {1}
        file_string=varargin{1};
        header_only=0;
        full_record=1;
    case {2}
        file_string=varargin{1};
        header_only=varargin{2}(1);
        full_record=1;
    case {3}
        file_string=varargin{1};
        header_only=varargin{2}(1);
        full_record=0;
        time_points_read=varargin{3}(1);
    otherwise
        error('Invalid number (%i) of arguments for load_slcontrol_file()',no_of_arguments);
end

% Variables

max_search_length=200;             % How far to look into the file before giving up
                                    % Saves scanning huge data files in their entirety if the
                                    % tag doesn't exist (it should be near
                                    % the beginning if it is in the header)

% Open data file
input_file=fopen(file_string,'r');
if (input_file == -1)
    error('Error: Data file %s could not be opened',file_string);
end
data_file.file_string=file_string;

% Skip the header if header_only = -1

if (header_only>=0)
    % Extract header file information using sub-functions below

    data_file.version=extract_string_value(input_file,max_search_length,'_VERSION_');
    data_file.file_info_string=extract_file_info(input_file,max_search_length);
    data_file.pCa=extract_float_from_string(data_file.file_info_string,'pCa');
    data_file.pH=extract_float_from_string(data_file.file_info_string,'pH');
    data_file.ADP=extract_float_from_string(data_file.file_info_string,'ADP');
    data_file.Pi=extract_float_from_string(data_file.file_info_string,'Pi');
    
    data_file.force_gain = extract_float_from_string(data_file.file_info_string,'Force_gain');
    
    data_file.drug = extract_string_value(input_file,max_search_length,'Drug');
    data_file.drug_concentration = extract_float_from_string(data_file.file_info_string,'Concentration_M');
    data_file.vehicle = extract_string_value(input_file,max_search_length,'Vehicle');
    data_file.vehicle_per_cent = extract_float_from_string(data_file.file_info_string,'Volume_per_cent');
    
    data_file.freeform_file = extract_string_value(input_file,max_search_length,'freeform_file:');

    data_file.creation_time.day=extract_integer_value(input_file,max_search_length,'Day');
    data_file.creation_time.month=extract_integer_value(input_file,max_search_length,'Month');
    data_file.creation_time.year=extract_integer_value(input_file,max_search_length,'Year');
    [data_file.creation_time.hours,data_file.creation_time.mins,data_file.creation_time.secs]= ...
        extract_time(input_file,max_search_length,'Time');

    data_file.FL_COMMAND=1.0e-6*extract_float_value(input_file,max_search_length,'FL_COMMAND_CALIBRATION:');
    data_file.FL_RESPONSE=1.0e-6*extract_float_value(input_file,max_search_length,'FL_RESPONSE_CALIBRATION:');
    data_file.FORCE_CALIBRATION=extract_float_value(input_file,max_search_length,'FORCE_CALIBRATION:');
    data_file.FL_POLARITY=extract_float_value(input_file,max_search_length,'FL_POLARITY:');
    data_file.SL_CALIBRATION=extract_float_value(input_file,max_search_length,'SL_volts_to_FL_COMMAND_volts:');

    data_file.muscle_length=1.0e-6*extract_float_value(input_file,max_search_length,'Muscle_length_(µm):');
    data_file.sarcomere_length=1.0e-9*extract_float_value(input_file,max_search_length,'Sarcomere_length_(nm):');
    data_file.area=extract_float_value(input_file,max_search_length,'Area_(m^2):');

    data_file.sampling_rate=extract_float_value(input_file,max_search_length,'sampling_rate_(Hz):');
    data_file.no_of_points=extract_integer_value(input_file,max_search_length,'data_points:');
    data_file.ktr_step=1.0e-6*extract_float_value(input_file,max_search_length,'ktr_step_(µm):');
    data_file.ktr_duration=1.0e-3*extract_float_value(input_file,max_search_length,'ktr_duration_(ms):');
    data_file.pre_ktr=1.0e-3*extract_float_value(input_file,max_search_length,'pre_ktr_(ms):');
    data_file.ktr_initiation_time=1.0e-3*extract_float_value(input_file,max_search_length,'ktr_initiation_time_ms:');
    data_file.post_ktr_servo_mode=extract_integer_value(input_file,max_search_length,'post_ktr_servo_mode:');

    data_file.no_of_triangles=extract_integer_value(input_file,max_search_length,'no_of_triangles:');
    data_file.ramp_mode=extract_integer_value(input_file,max_search_length,'ramp_mode:');
    data_file.hold_mode=extract_integer_value(input_file,max_search_length,'hold_mode:');
    data_file.hold_time=1e-3*extract_float_value(input_file,max_search_length,'hold_(ms):');
    data_file.triangle_size=1.0e-6*extract_float_value(input_file,max_search_length,'triangle_size_(µm):');
    data_file.pre_triangle_time=1.0e-3*extract_float_value(input_file,max_search_length,'pre_triangle_(ms):');
    data_file.triangle_halftime=1.0e-3*extract_float_value(input_file,max_search_length,'triangle_halftime_(ms):');
    data_file.inter_triangle_interval=1.0e-3*extract_float_value(input_file,max_search_length,'inter_triangle_interval_(ms):');
    data_file.relative_first_triangle_size=extract_float_value(input_file,max_search_length,'relative_first_tri_size:');

    data_file.first_step_time=1.0e-3*extract_float_value(input_file,max_search_length,'first_step_ms:');
    data_file.first_step_size=1.0e-6*extract_float_value(input_file,max_search_length,'first_step_microns:');
    data_file.second_step_time=1.0e-3*extract_float_value(input_file,max_search_length,'second_step_ms:');
    data_file.second_step_size=1.0e-6*extract_float_value(input_file,max_search_length,'second_step_microns:');
    
    % FV parameters
    data_file.pre_servo_time = 1.0e-3*extract_float_value(input_file,max_search_length,'pre_servo_(ms):');
    data_file.servo_time = 1.0e-3*extract_float_value(input_file,max_search_length,'servo_(ms):');
    data_file.isotonic_hold = extract_float_value(input_file,max_search_length,'isotonic_hold_(V):');
 
    
    data_file.no_of_channels=extract_integer_value(input_file,max_search_length,'No_of_input_channels:');
    if (isnan(data_file.no_of_channels))
        data_file.no_of_channels=4;
    end
end

% Now read channel data if header_only~=1
% Scan the file until you find the data tag

if header_only~=1
	frewind(input_file);
	counter=1;
	test_string=fscanf(input_file,'%s',1);
	while ((strcmp(test_string,'Data')~=1)&&(counter<=max_search_length))
        test_string=fscanf(input_file,'%s',1);
        counter=counter+1;
	end
	
	if (counter>max_search_length)
        error('No "Data" tag found in file %s\nIs this an SLControl data file?',file_string);
        return
	end
	
	% Now read the data and store as structure elements
	
    if (full_record)
    	data=fscanf(input_file,'%f',[data_file.no_of_channels+1,inf]);
    else
        data=fscanf(input_file,'%f', ...
            [data_file.no_of_channels+1,time_points_read]);
    end
	
    if (data_file.no_of_channels==4)
        data_file.time=(data(1,:))';
        data_file.force=(data(2,:))';
        data_file.sl=(data(3,:))';
        data_file.fl=(data(4,:))';
        data_file.intensity=(data(5,:))';
    end
    
    if (data_file.no_of_channels==8)
        data_file.time=(data(1,:))';
        data_file.force=(data(2,:))';
        data_file.fl=(data(3,:))';
        data_file.fl_command=(data(4,:))';
        data_file.stimulus=(data(5,:))';
        data_file.light_intensity=(data(6,:))';
        data_file.shutter=(data(7,:))';
        data_file.wavelength_command=(data(8,:))';
        data_file.spare=(data(9,:))';
    end
    
    % Adjust for force_gain
    if (~isnan(data_file.force_gain))
        data_file.force = data_file.force ./ data_file.force_gain;
    end
    
end

% Close the input_file
fclose(input_file);

% Sub-functions

% Extract integer after tag_string

function integer=extract_integer_value(input_file,max_search_length,tag_string);

frewind(input_file);
counter=1;
test_string=fscanf(input_file,'%s',1);
while ((strcmp(test_string,tag_string)~=1)&&(counter<=max_search_length))
    test_string=fscanf(input_file,'%s',1);
    counter=counter+1;
end
if (counter<=max_search_length)
    integer=fscanf(input_file,'%i',1);
else
    integer=NaN;
end


% Extract float after tag_string

function float=extract_float_value(input_file,max_search_length,tag_string);

frewind(input_file);
counter=1;
test_string=fscanf(input_file,'%s',1);
while ((strcmp(test_string,tag_string)~=1)&&(counter<=max_search_length))
    test_string=fscanf(input_file,'%s',1);
    counter=counter+1;
end
if (counter<=max_search_length)
    float=fscanf(input_file,'%f',1);
else
    float=NaN;
end

% Extract string after tag_string

function string=extract_string_value(input_file,max_search_length,tag_string);

frewind(input_file);
counter=1;
test_string=fscanf(input_file,'%s',1);
while ((strcmp(test_string,tag_string)~=1)&&(counter<=max_search_length))
    test_string=fscanf(input_file,'%s',1);
    counter=counter+1;
end
if (counter<=max_search_length)
    string=fscanf(input_file,'%s',1);
else
    string='Default value';
end

% Extract all strings between _FILE_INFO_START_ and _FILE_INFO_STOP_ and
% return as a single string

function string=extract_file_info(input_file,max_search_length);

string='';
frewind(input_file);
counter=1;
test_string=fscanf(input_file,'%s',1);
while ((strcmp(test_string,'_FILE_INFO_START_')~=1)&&(counter<=max_search_length))
    test_string=fscanf(input_file,'%s',1);
    counter=counter+1;
end
if (counter<=max_search_length)
    counter=1;
    temp_string='';
    while (strcmp(temp_string,'_FILE_INFO_STOP_')~=1)
        temp_string=fscanf(input_file,'%s',1);
        if (strcmp(temp_string,'_FILE_INFO_STOP_')~=1)&(counter>1)
            string=sprintf('%s %s',string,temp_string);
        else
            if (strcmp(temp_string,'_FILE_INFO_STOP_')~=1)
                string=sprintf('%s',string,temp_string);
            end
        end
        counter=counter+1;
    end
else
    string='';
end

% Extract time after tag_string

function [hr,min,sec]=extract_time(input_file,max_search_length,tag_string);

frewind(input_file);
counter=1;
test_string=fscanf(input_file,'%s',1);
while ((strcmp(test_string,tag_string)~=1)&&(counter<=max_search_length))
    test_string=fscanf(input_file,'%s',1);
    counter=counter+1;
end
if (counter<=max_search_length)
    hr=fscanf(input_file,'%i',1);
    test_string=fscanf(input_file,'%s',1);
    min=fscanf(input_file,'%i',1);
    test_string=fscanf(input_file,'%s',1);
    sec=fscanf(input_file,'%i',1);
else
    hr=NaN;
    min=NaN;
    sec=NaN;
end

% Extract values from file_info_string

function number=extract_float_from_string(file_info_string,tag_string);

start_search=findstr(file_info_string,tag_string)+length(tag_string);
if (start_search)
    number=str2num(char(strtok(file_info_string(start_search:length(file_info_string)))));
    if (isempty(number))
        number = NaN;
    end
else
    number=NaN;
end


