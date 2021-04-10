function bootstrap_dec_rval_corr(D,nReps,do_trim,do_shuffle)
% bootstrap the distribution of r-values computed on the correlation of 
% the correlation coefficients between to matrices

pairs = unique_idx_pairs(numel(D),0);

% initialize r-value placeholder
r_bs = NaN(nReps,size(pairs,1));

% iterate over replicates and pairs of data sets
fprintf('\n');
for i=1:nReps
    fprintf('Iteration %i of %i\n',i,nReps)
    for j=1:size(pairs,1)
        % Bootstrap the data
        DA_bs = D(pairs(j,1));
        da = DA_bs.data;
        D1_idx = randi(size(da,1),[size(da,1) 1]);
        DA_bs.data = da(D1_idx,:);
        DB_bs = D(pairs(j,2));
        db = DB_bs.data;
        D2_idx = randi(size(db,1),[size(db,1) 1]);
        DB_bs.data = db(D2_idx,:);
        [DA_bs, DB_bs] = pairFields(DA_bs,DB_bs,'Trim',do_trim);
        
        if do_shuffle
            DA_bs.data = shuffle_columns(DA_bs.data);
            DB_bs.data = shuffle_columns(DB_bs.data);
        end
        [r,~] = corr_of_corrcoef(DA_bs, DB_bs, 'Plot', false);
        r_bs(i,j) = r(1,2);
    end
end


ah = gca;
line_handles = [];
for i=1:size(pairs,1)
    lh = plot_kde(r_bs(:,i),ah,[-1 1],ah.ColorOrder(i,:));
    line_handles = [line_handles; lh];
end

labels = cellfun(@(p) sprintf('D%i-D%i',p(1),p(2)), num2cell(pairs,2),...
    'UniformOutput', false);
legend(line_handles,labels,'Location','Northwest');
xlabel('correlation coefficient');

