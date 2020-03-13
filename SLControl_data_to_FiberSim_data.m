function SLControl_data_to_MyoSim_data(varargin)

% Take inputs
p=inputParser;
addParamValue(p,'slc_data',[]);
addParamValue(p,'start_time_s',0);
addParamValue(p,'stop_time_s',[]);
addParamValue(p,'t_inc',[]);
addParamValue(p,'t_zoom',[]);
addParamValue(p,'pCa_data',[]);
addParamValue(p,'half_sarcomere_length_nm',[]);
addParamValue(p,'pre_data',[]);
addParamValue(p,'protocol_file_string',[]);
addParamValue(p,'target_file_string',[]);
addParamValue(p,'smooth_fl',0);
addParamValue(p,'slc_figure_number',0);
addParamValue(p,'template_figure_number',0);

parse(p,varargin{:});

field_names = fieldnames(p.Results);
for i=1:numel(field_names)
    eval_string = sprintf('%s = p.Results.%s;', ...
        field_names{i},field_names{i});
    eval(eval_string);
end

% Code

td = p.Results.slc_data;

% Some updates
if (isempty(stop_time_s))
    % Load first file
    stop_time_s = d.time(end);
end

% Create holders
dt = [];
target_force = [];
holder_dhsl = [];
holder_pCa = [];

% Pull data

% Interpolate to different dt if required
t = td.time;
if (~isempty(t_inc))
    if (isempty(t_zoom))
        td_new.time = 0:t_inc:t(end);
    else
        td_new.time = [(0:t_inc:(t_zoom(1)-t_inc))' ;
            (t_zoom(1):(0.1*t_inc):t_zoom(2))' ;
            ((t_zoom(2)+t_inc):t_inc:t(end))']';
    end
    
    td_new.force = interp1(td.time,td.force,td_new.time);
    td_new.fl = interp1(td.time,td.fl,td_new.time);

    td.time = td_new.time';
    td.force = td_new.force';
    td.fl = td_new.fl';
end

ti = find(td.time > start_time_s,1,'first') : ...
        find(td.time < stop_time_s,1,'last');
    
ddt = diff(td.time(ti));
ddt = [ddt ; ddt(end)];
        
% Force
target_force = [target_force td.force(ti)];

% FL
nfl = td.fl(ti)./td.fl(ti(i));
if (smooth_fl)
    nfl2 = movmean(nfl,smooth_fl)
end

figure(3);
clf
hold on;
plot(nfl,'b-');
plot(nfl2,'r-');


nfl = nfl2;
dnfl = diff(nfl);
dnfl = [dnfl ; dnfl(end)];
dhsl = dnfl * half_sarcomere_length_nm;

holder_dhsl = [holder_dhsl dhsl];
    
% pCa
if (isempty(pCa_data))
    holder_pCa = [holder_pCa td.pCa*ones(numel(ti),1)];
end

% Add in pre_data
if (~isempty(pre_data))
    sa = size(pre_data.dt)
    sb = size(ones(numel(ti),1))
    dt = [pre_data.dt ; ddt];
    holder_dhsl = [pre_data.dhsl ; holder_dhsl];
    holder_pCa = [pre_data.pCa ; holder_pCa];
end

% Output
sdt = size(dt)
sp = size(holder_pCa)
sdhsl = size(holder_dhsl)


holder_mode = -2 * ones(numel(dt),1);
out_t = table(dt, holder_pCa, holder_dhsl, holder_mode);
writetable(out_t, protocol_file_string,'Delimiter','\t');

out_f = table(target_force);
writetable(out_f, target_file_string, 'Delimiter','\t');

% Draw template
if (template_figure_number)
    
    t = cumsum(dt);
    
    figure(template_figure_number);
    clf;
    
    subplot(4,1,1);
    hold on;
    [r,c]=size(target_force);
    if (isempty(pre_data))
        vi = 1:numel(dt);
    else
        vi = (numel(pre_data.dt(:,1))+1):numel(dt);
    end
    
    for i=1:c
        plot(t(vi),target_force(:,i));
    end
    xlim([0 t(end)]);
    
    subplot(4,1,2);
    hold on;
    [r,c]=size(holder_dhsl);
    for i=1:c
        plot(t,cumsum(holder_dhsl(:,i)));
    end
    xlim([0 t(end)]);

    subplot(4,1,3);
    hold on;
    [r,c]=size(holder_pCa);
    for i=1:c
        plot(t,holder_pCa(:,i));
    end
    xlim([0 t(end)]);

    subplot(4,1,4);
    hold on;
    dd = diff(holder_dhsl);
    dd = [dd ; dd(end)];

    for i=1:c
        plot(t,dd);
    end
    xlim([0 t(end)]);

end






