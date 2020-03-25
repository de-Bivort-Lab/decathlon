

for i=1:numel(D)
    % find clumpiness measures
    [fields, p] = groupFields(standardize_fieldnames(D(i).fields), 'clumpiness');
    
    % exclude bout clumpiness measures
    [assay,metric] = parse_fieldnames(fields);
    keep_fields = ~str_list_contains(metric,'bout');
    clump_idx = p(keep_fields);
    clump_fields = cellfun(@(a,m) sprintf('%s %s',a,m), ...
        assay(keep_fields), metric(keep_fields),'UniformOutput', false);
    [clump_fields,ii] = sort(clump_fields);
    clump_idx = clump_idx(ii);
    clump_assay = parse_fieldnames(clump_fields);
    
    % find clumpiness measures
    [fields, p] = groupFields(standardize_fieldnames(D(i).fields), 'switchiness');
    
    % exclude bout clumpiness measures
    [assay,metric] = parse_fieldnames(fields);
    keep_fields = ~str_list_contains(metric,'bout');
    switch_idx = p(keep_fields);
    switch_fields = cellfun(@(a,m) sprintf('%s %s',a,m), ...
        assay(keep_fields), metric(keep_fields),'UniformOutput', false);
    [switch_fields,ii] = sort(switch_fields);
    switch_idx = switch_idx(ii);
    switch_assay = parse_fieldnames(switch_fields);
    
    assay_matches = str_list_contains(clump_assay,switch_assay);
    keep_clump = any(assay_matches,2);
    keep_switch = any(assay_matches);
    switch_idx = switch_idx(keep_switch);
    switch_fields = switch_fields(keep_switch);
    switch_assay = switch_assay(keep_switch);
    clump_idx = clump_idx(keep_clump);
    clump_fields = clump_fields(keep_clump);
    clump_assay = clump_assay(keep_clump);
    
    % plot pairs
    nplots = max([numel(clump_idx) numel(switch_idx)]);
    pair_idx = cellfun(@(sa) find(strcmpi(clump_assay,sa),1), switch_assay);
    pair_idx = [clump_idx(pair_idx)' switch_idx'];
    scatter_field_pairs(D(i),pair_idx)
end