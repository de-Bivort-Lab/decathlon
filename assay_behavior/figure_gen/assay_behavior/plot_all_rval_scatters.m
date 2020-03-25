

% standardize batches
D = standardize_data_by_group(D);

% standardize metrics
D = standardize_by_field(D);

% impute the data
D = impute_decathlon_structs(D,'ImputeMode','mean');

% define dataset pairs
pairs = unique_idx_pairs(numel(D),1);

figure;
for i=1:4
    switch i
        case 1
            collapse_mode = 'average';
            collapse_fields = 'none';
            title_str = 'circ expanded';
        case 2
            collapse_mode = 'average';
            collapse_fields = 'circadian';
            title_str = 'circ collapsed';
        case 3
            collapse_mode = 'average';
            collapse_fields = 'all';
            title_str = 'a priori collapsed';
        case 4
            collapse_mode = 'PCA';
            collapse_fields = 'all';
            title_str = 'a priori PCA';
    end
    opts = {'CollapseFields';collapse_fields;'CollapseMode';collapse_mode;...
        'PCs';2;'Trim';true};
    D_col = pair_decathlon_structs(D,opts{:});
    for j=1:size(pairs,1)
        subplot(4,numel(D_col),(i-1)*size(pairs,1)+j);
        plot_title = {'';sprintf('D%i-D%i (%s)',pairs(j,1),pairs(j,2),title_str)};
        corr_of_corrcoef(D_col(pairs(j,1)),D_col(pairs(j,2)),...
            'Title',plot_title,'Plot',true);
        axis equal
    end
end

