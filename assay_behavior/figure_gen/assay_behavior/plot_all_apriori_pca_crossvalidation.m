
% set file path to decathlon data
fdir = '';
D = load_decathlon_structs(fdir,'D13_als');
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
        min_idx{i}(j) = 1;
         if size(apriori_data{j},2)>1
            [train,test] = cross_validate_pca(apriori_data{j},'KFolds',50,'TestSize',0.25);
            [~,min_idx{i}(j)] = min(mean(test,2));
        end
        title(sprintf('D%i - %s, PC thresh = %i',i,apriori_names{j},min_idx{i}(j)));
    end
end