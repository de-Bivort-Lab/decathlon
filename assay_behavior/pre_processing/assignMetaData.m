function meta = assignMetaData(labels_table, name, date)


plate_map = cell(192*2,1);
plate_map(1:96) = {'A'};
plate_map(97:192) = {'B'};
plate_map(193:288) = {'C'};
plate_map(289:384) = {'D'};
num_traces = sum(~isnan(labels_table.ID));
meta.Day = repmat(labels_table.Day(1),num_traces,1);

% calculate time of day from date
date(date=='_')=[];
b = find(date=='-');
tod = str2double(date(b(3)+1:b(4)-1))*3600 +...
    str2double(date(b(4)+1:b(5)-1))*60 + str2double(date(b(5)+1:end));

% assign time of day and plate meta data
meta.TimeofDay = repmat(tod,num_traces,1);

switch name
    case 'Olfaction'
        meta.Plate = plate_map(labels_table.ID);
        
    otherwise

        meta.Plate = plate_map(labels_table.ID);

        switch name
            case 'Circadian'
                meta.Box = arrayfun(@(x) ['circ-' num2str(x)],...
                    labels_table.Box,'UniformOutput',false);

            otherwise
                meta.Box = arrayfun(@(x) ['proj-' num2str(x)],...
                    labels_table.Box,'UniformOutput',false);
                meta.Tray = num2cell(num2str(labels_table.Tray),2);
        end
end
        
        