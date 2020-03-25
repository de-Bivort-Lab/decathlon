function fields = standardize_fieldnames(fields)

% convert Day 1 Y-maze entries to Culling
[grp_f, fidx] = groupFields(fields, '(1)');
[assays, metrics, days] = parse_fieldnames(grp_f);
filt = strcmpi(assays,'Y-Maze');
new_fields = arrayfun(@(m,d) sprintf('Culling %s (%i)',m{1},d),...
    metrics(filt), days(filt), 'UniformOutput', false);
new_fields = remove_NaN_days(new_fields, isnan(days(filt)));
fields(fidx(filt)) = new_fields;

% convert Arena Circling to Arena
[assays, metrics, days] = parse_fieldnames(fields);
filt = strcmpi(assays,'Arena Circling') | strcmpi(assays,'Basic Tracking');
new_fields = arrayfun(@(m,d) sprintf('Arena %s (%i)',m{1},d),...
    metrics(filt), days(filt), 'UniformOutput', false);
new_fields = remove_NaN_days(new_fields, isnan(days(filt)));
fields(filt) = new_fields;

% convert circling_mu to circling
[assays, metrics, days] = parse_fieldnames(fields);
filt = strcmpi(metrics,'circling_mu');
new_fields = arrayfun(@(m,d) sprintf('%s circling (%i)',m{1},d),...
    assays(filt), days(filt), 'UniformOutput', false);
new_fields = remove_NaN_days(new_fields, isnan(days(filt)));
fields(filt) = new_fields;

% convert Y-maze to Y-Maze
[assays, metrics, days] = parse_fieldnames(fields);
filt = strcmpi(assays,'Y-Maze');
new_fields = arrayfun(@(m,d) sprintf('Y-Maze %s (%i)',m{1},d),...
    metrics(filt), days(filt), 'UniformOutput', false);
new_fields = remove_NaN_days(new_fields, isnan(days(filt)));
fields(filt) = new_fields;

% convert Y-maze to Y-Maze
[assays, metrics, days] = parse_fieldnames(fields);
filt = strcmpi(assays,'LED Y-Maze');
new_fields = arrayfun(@(m,d) sprintf('LED Y-Maze %s (%i)',m{1},d),...
    metrics(filt), days(filt), 'UniformOutput', false);
new_fields = remove_NaN_days(new_fields, isnan(days(filt)));
fields(filt) = new_fields;

% convert pre odor to preodor
filt = strcmpi(metrics,'preodor_occupancy') | strcmpi(metrics,'pre_odor_occupancy');
new_fields = arrayfun(@(d) sprintf('Olfaction %s (%i)','preodor_occupancy',d),...
    days(filt), 'UniformOutput', false);
new_fields = remove_NaN_days(new_fields, isnan(days(filt)));
fields(filt) = new_fields;


function new_fields = remove_NaN_days(new_fields, nan_days)

day_label_removed = regexp(new_fields(nan_days),' \(.*\)','split');
day_label_removed = cellfun(@(dlr) dlr{1}, day_label_removed, 'UniformOutput', false);
new_fields(nan_days) = day_label_removed;