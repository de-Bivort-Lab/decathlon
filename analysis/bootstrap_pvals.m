function pval_bs = bootstrap_pvals(D,nReps)

% initialize r-value placeholder
pair_idx = upper_triangle_idx(numel(D.fields));
pval_bs = NaN(numel(pair_idx),nReps);


% iterate over replicates and pairs of data sets
fprintf('\n');
for i=1:nReps
    
    fprintf('Iteration %i of %i\n',i,nReps)

    % randomly sample data with replacement
    bs_idx = randi(size(D.data,1),[size(D.data,1) 1]);
    data = D.data(bs_idx,:);
    
    % compute correlation pvals and store pvals for unique pairs
    [~,p_tmp] = corr(data,'Type','spearman','rows','pairwise');
    pval_bs(:,i) = p_tmp(pair_idx);
end


plot_cdf_patch(pval_bs(:));

% compute correlation pvals and store pvals for unique pairs
[~,p_tmp] = corr(D.data,'Type','spearman','rows','pairwise');
plot_cdf_patch(p_tmp(pair_idx));


    