
D = D13_als;
max_reps = 200;
d_labels = {'D12';'D3'};

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
    
    % create pairs data structs with current options
    opts = {'CollapseFields';collapse_fields;'CollapseMode';collapse_mode;...
        'PCs';2;'Trim';true;'ImputeMode';'none'};
    D_p = D;
    D_p = pair_decathlon_structs(D_p,opts{:});
    D_p = standardize_by_field(D_p);
    
    for j=1:numel(D)
        % open new subplot
        subplot(4,numel(D),(i-1)*numel(D)+j);
        plot_pca_n_drop_scree(D_p(j).data,max_reps);
        title(sprintf('%s (%s)',d_labels{j},title_str));
    end
end

