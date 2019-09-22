function [eigvecs, eigvals, var_explained] = bootstrap_pca_nullmodel(data,nreps,max_n)

% check that matrix is complete
if any(isnan(data(:)))
    error('Cannot perform PCA. Matrix contains missing values.');
end

eigvecs = cell(nreps,1);
eigvals = cell(nreps,1);
var_explained = cell(nreps,1);
for i=1:nreps
    resampled_data = data(randi(size(data,1),[max_n 1]),:);
    shuffled_data = shuffle_columns(resampled_data);
    [eigvecs{i}, ~, eigvals{i}, ~, var_explained{i}] = pca(shuffled_data);
end