function [assays, metrics, days] = parse_fieldnames(fields)
% Parse decathlon behavior field names from 'Assay metric (day)' format 
% into component parts

% get assay name
names = regexp(fields,'(\<[A-Z][\w|-]*)*','match');
assays = cell(numel(names),1);
for i=1:numel(assays)
    assays{i} = names{i}{1};
    for j=2:numel(names{i})
        assays{i} = sprintf('%s %s',assays{i},names{i}{j});
    end
end

% get metric name
metrics = regexp(fields,'(?<!-)(\<[a-z][\w|_]*)*','match');
metrics = cat(1,metrics{:});

% get day number
days = regexp(fields,'(?<=\()[0-9]*(?=\))','match');
days(cellfun(@isempty,days)) = {{''}};
days = cellfun(@(d) d{1}, days, 'UniformOutput', false);
days = cellfun(@str2double, days);