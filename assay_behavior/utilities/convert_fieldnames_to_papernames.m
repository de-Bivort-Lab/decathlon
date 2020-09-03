function new_fields = convert_fieldnames_to_papernames(fields)

[assay,metric,day] = parse_fieldnames(fields);

% assay name conversion list
input_assay= {'culling';'slow phototaxis';'temporal phototaxis';'Olfaction'};
output_assay = {'Activity Screen';'Spatial Shade-light';'Temporal Shade-light';'Odor Sensitivity'};
for i=1:numel(input_assay)
   assay(strcmpi(assay,input_assay{i})) = output_assay(i); 
end

input_metric = {'circling';'nTrials';'nBouts';'right_bias';...
    'light_bias';'hand_clumpiness';'hand_switchiness';'speed';...
    'occupancy'};
output_metric = {'circling_bias';'choice_number';'bout_number';'turn_bias';...
    'photo_bias';'turn_clumpiness';'turn_switchiness';'mean_speed';...
    'photo_occupancy'};
for i=1:numel(input_metric)
   metric(strcmpi(metric,input_metric{i})) = output_metric(i); 
end

new_fields = cellfun(@(a,m) sprintf('%s %s',a,m), assay, metric, 'UniformOutput', false);

% handle special cases
input_field = {'Optomotor choice_number';'Spatial Shade-light choice_number';...
    'Temporal Shade-light choice_number';'Odor Sensitivity occupancy'};
output_field = {'Optomotor trial_number';'Spatial Shade-light trial_number';...
    'Temporal Shade-light trial_number';'Odor Sensitivity odor_occupancy'};
for i=1:numel(input_field)
   new_fields(strcmpi(new_fields,input_field{i})) = output_field(i); 
end

for i=1:numel(fields)
    if any(fields{i}=='(')
        new_fields{i} = sprintf('%s (%i)',new_fields{i},day(i));
    end
end