function norm_data = quantile_normalize(data)
% perform quantile normalization on each column

% sort each column in ascending order
[sorted_cols,row_perm] = cellfun(@sort,num2cell(data,1),'UniformOutput', false);
sorted_cols = cat(2,sorted_cols{:});
%[~,~,row_ranks] = cellfun(@unique,num2cell(data,1),'UniformOutput', false);

row_ranks = NaN(size(data));
for i=1:numel(row_perm)
   row_ranks(row_perm{i},i) = 1:numel(row_perm{i});
end

% calculate row-wise means
row_means = nanmean(sorted_cols,2);

norm_data = NaN(size(data));
for i=1:size(row_ranks,2)
    norm_data(:,i) = row_means(row_ranks(:,i));
end
norm_data(isnan(data)) = NaN;

