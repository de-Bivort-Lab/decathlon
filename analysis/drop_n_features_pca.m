function var_exp = drop_n_features_pca(data,n,max_reps,shuffle)
% iteratively drop features from data and perform pca

% generate random index groupings
max_n = size(data,2);
idx_grps = arrayfun(@(ii) randperm(max_n,n), 1:max_reps*3, 'UniformOutput', false);
unique_idx_grps = cell(max_reps,1);

% iterate over random groups and select unique groups
i = 1;
j = 1;
while j<=numel(unique_idx_grps) && i < numel(idx_grps)
    tmp_grp = idx_grps{i};
    dup_grp = any(cellfun(@(uig) all(ismember(tmp_grp,uig)), unique_idx_grps));
    
    % if group is not duplicate, add it to grp list
    if ~dup_grp
        unique_idx_grps{j} = tmp_grp;
        j=j+1;
    end
    i=i+1;
end

% remove empty grps
unique_idx_grps(cellfun(@isempty,unique_idx_grps))=[];

% do PCA, iteratively dropping out groups
var_exp = NaN(numel(unique_idx_grps),max_n-n);
for i=1:numel(unique_idx_grps)
    filt_data = data;
    filt_data(:,unique_idx_grps{i}) = [];
    if shuffle
        filt_data = shuffle_columns(filt_data);
    end
    [~,~,~,~,var_exp(i,:)] = pca(filt_data,'NumComponents',size(filt_data,2));
end

