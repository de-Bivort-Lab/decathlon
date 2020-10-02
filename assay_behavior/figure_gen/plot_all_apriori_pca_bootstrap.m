
% set file path to decathlon data
fdir = pwd;
D = load_decathlon_structs(fdir,'D_als_filled_batch_merged');
D_p = pair_decathlon_structs(D);

% cross-validate PCA for all apriori groups
nrow = 3;
ncol = 4;
min_idx = cell(numel(D_p),1);
for i=1:numel(D_p)
    figure;
    [apriori_data, apriori_names] = group_apriori_fields(D_p(i));
    min_idx{i} = NaN(numel(apriori_data),1);
    for j=1:numel(apriori_data)
        subplot(nrow,ncol,j);
        if size(apriori_data{j},2)>1
            [~,~,nkeep]=plot_pca_bootstrap(apriori_data{j},50,95,'noncummulative');
        end
        title(sprintf('%s, k=%i',apriori_names{j},nkeep));
    end
end