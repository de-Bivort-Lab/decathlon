%% Get parent directory of all decathlon files

root = 'D:\decathlon_data\decathlon_data\';
sr = 'decathlon 2-2018\D1D2_data\';
fdir = {sprintf('%s%s%s',root,sr,'D1_data\');...
    sprintf('%s%s%s',root,sr,'D2_data\');...
    sprintf('%s%s',root,'decathlon 1-2019\data\')};
load('D:\decathlon_data\decathlon_data\decathlon 1-2019\meta\culling_ID_permutation.mat');

%% build the decathlon data struction from individual assay files

fprintf('\n');
for i=1:numel(fdir)
    fprintf('Building struct %i of %i\n',i,numel(fdir))
    if i==3
        options = {'FieldFunction',@getDataFields,'CullingPermutation',ID};
    else
        options = {'FieldFunction',@getDataFields_legacy,'CullingPermutation',[]};
    end
    D(i) = build_DecathlonDataStruct(fdir{i},options{:});
end

%% pre-process and format data matrices

clearvars -except D
D_raw = D;
D12_cat = [cat_decathlon_structs(D_raw(1:2),...
    'ImputeMode','none','Standardize',false); D_raw(3)];
D_z = impute_decathlon_structs(D,'ImputeMode','none');
[D123_als, als_data] = avg_als_impute(D_z,200);
D123_als = standardize_by_field(D123_als);
D12_cat = cat_decathlon_structs(D123_als(1:2),'ImputeMode','none','Standardize',false);
D12_cat_als = avg_als_impute(D12_cat,200);
D12_cat_als = standardize_by_field(D12_cat_als);
D13_als = [D12_cat_als D123_als(3)];
