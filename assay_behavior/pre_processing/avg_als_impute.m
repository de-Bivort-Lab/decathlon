function [D, als_data] = avg_als_impute(D, nReps)

% initialize placeholder for data
als_data = cell(numel(D),1);
for i=1:numel(D)
    als_data{i} = NaN([size(D(i).data) nReps]);
end

% generate nReps als imputed data sets for each decathlon struct
fprintf('\n');
for i=1:nReps
    fprintf('ALS imputation: iteration %i of %i\n',i,nReps);
    D_imputed = impute_decathlon_structs(D,'ImputeMode','als','Standardize',false);
    for j=1:numel(D)
        als_data{j}(:,:,i) = D_imputed(j).data;
    end
end

% average all runs together
als_avg = cellfun(@(d) mean(d,3), als_data, 'UniformOutput',false);
als_med= cellfun(@(d) median(d,3), als_data, 'UniformOutput',false);
for i=1:numel(D)
    D(i).imputed = isnan(D(i).data);
    D(i).data = als_med{i};  
end