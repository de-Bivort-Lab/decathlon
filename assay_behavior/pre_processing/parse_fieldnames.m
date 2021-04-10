function [assays, metrics, days] = parse_fieldnames(fields,varargin)
% Parse decathlon behavior field names from 'Assay metric (day)' format 
% into component parts

trim = false;
if ~isempty(varargin)
    trim = varargin{1};
end
assay_long= {'Circadian';'Culling';'LED Y-Maze';'Olfaction';'Optomotor';...
    'Slow Phototaxis';'Temporal Phototaxis'};
assay_short= {'Circ.';'Cull.';'LED-Y.';'Olf.';'Opto.';'Spat. SL';'Temp. SL'};

metric_long = {'bout_length';'bout_clumpiness';'circling_blank';'occupancy';...
    'optomotor_index';'right_bias';'preodor_occupancy';'hand_clumpiness';...
    'hand_switchiness';'light_bias';'light_switchiness';...
    'gravitactic_index';'light_clumpiness'};
metric_short = {'b_len';'b_clump';'circ_blank';'occ';'opto_ind';'r_bias';...
    'preodor_occ';'hand_clump';'hand_switch';'lit_bias';'lit_switch';'grav_ind';'lit_clump'};

% get assay name
names = regexp(fields,'(\<[A-Z][\w|-]*)*','match');
assays = cell(numel(names),1);
for i=1:numel(assays)
    assays{i} = names{i}{1};
    for j=2:numel(names{i})
        assays{i} = sprintf('%s %s',assays{i},names{i}{j});
    end
end
if trim
    for i=1:numel(assay_long)
       assays(strcmp(assays,assay_long{i})) = assay_short(i); 
    end
end

% get metric name
metrics = regexp(fields,'(?<!-)(\<[a-z][\w|_]*)*','match');
metrics = cat(1,metrics{:});
if trim
    for i=1:numel(metric_long)
       metrics(strcmp(metrics,metric_long{i})) = metric_short(i); 
    end
end

% get day number
days = regexp(fields,'(?<=\()[0-9]*(?=\))','match');
days(cellfun(@isempty,days)) = {{''}};
days = cellfun(@(d) d{1}, days, 'UniformOutput', false);
days = cellfun(@str2double, days);