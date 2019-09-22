function [apriori_grouped, apriori_group_names, grp_idx] = group_apriori_fields(D)


% define keywords for clustering groups
group_keywords = [{{'speed';'nTrials';'nBouts';'bout_length'}};...
    {{'bout_clumpiness'}};...
    {'Circadian gravi'};...
    {{'hand_clumpiness';'iti';'light_clumpiness'}};...
    {'switchiness'};...
    {{'Arena circling';'Phototaxis circling';'Culling circling';'right_bias'}};...
    {'Circadian circling'};...
    {{'light_bias';'Phototaxis occupancy'}};...
    {'Optomotor circling';'optomotor_index';'Olfaction occupancy'}];

% initialize grouped field names
apriori_group_names = {'Activity';'Bout clumpiness';'Gravitaxis';'Clumpiness';'Switchiness';'Handedness';...
                        'Circadian circling';'Phototaxis';...
                        'Opto Handedness';'Opto index';'Odor sensitivity'};

% get group indices and group data by indices
std_fn = standardize_fieldnames(D.fields);
[~, grp_idx] = arrayfun(@(gkw) groupFields(std_fn, gkw), group_keywords,...
                    'UniformOutput', false);
for i=1:numel(grp_idx)
    
   ii = grp_idx{i};
   tmp_fields = D.fields(ii);
   [assay,metric,day] = parse_fieldnames(tmp_fields);
   
   % sort first by assay
   unique_assays = sort(unique(assay));
   unique_assays(strcmpi(unique_assays,'circadian')) = [];
   unique_assays = ['Circadian';unique_assays];
   sorted_ii = NaN(numel(ii),1);
   for j=1:numel(unique_assays)
       
       jj = find(strcmpi(assay,unique_assays{j}));
       
       % unique metrics
       tmp_metrics = metric(jj);
       unique_metrics = sort(unique(tmp_metrics));
       sorted_jj = NaN(size(jj));
       for k=1:numel(unique_metrics)
           kk = find(strcmpi(tmp_metrics,unique_metrics{k}));
           j_ct = find(isnan(sorted_jj),1);
           sorted_jj(j_ct:j_ct+numel(kk)-1) = jj(kk);
       end
       
       i_ct = find(isnan(sorted_ii),1);
       sorted_ii(i_ct:i_ct+numel(jj)-1) = ii(sorted_jj);
   end
   
   grp_idx{i} = sorted_ii;
end
                
apriori_grouped = cellfun(@(idx) D.data(:,idx), grp_idx, 'UniformOutput', false);


