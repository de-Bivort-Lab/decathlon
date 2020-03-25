function [permutation, sorted_labels] = get_cluster_permutation(D)
% D = decathlon struct or array of decathlon structs

data = cat(1,D.data);

% calculate covariance matrix
[r,p] = corrcoef(data,'rows','pairwise');

% replace NaNs
r(isnan(r))=0;
p(isnan(p))=1;

% get linkage values
Z=linkage(r,'single','spearman');
f=figure;
[~, ~, permutation]=dendrogram(Z,length(r));
close(f);

% labels
labels = pretty_labels(D(1).fields);
sorted_labels = labels(permutation);