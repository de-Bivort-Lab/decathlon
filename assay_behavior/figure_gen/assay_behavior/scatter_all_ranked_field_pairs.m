% sort metric pairs by ranked p-value and scatter resulting top pairs

% define collapsing parameters
collapse_mode = 'PCA';
collapse_fields = 'all';
pcs = 2;
do_trim = true;
impute_mode = 'none';
opts = {'CollapseFields';collapse_fields;'CollapseMode';collapse_mode;...
    'PCs';pcs;'Trim';do_trim;'ImputeMode';impute_mode};

% standardize batches
D = standardize_data_by_group(D);

% standardize metrics
D = standardize_by_field(D);

% impute the data
D = impute_decathlon_structs(D,opts{:});

% plot correlation matrices
[fh1, fh2] = plot_paired_corrmats(D,opts{:});

% pair the struct fields
D = pair_decathlon_structs(D,opts{:});
D = standardize_by_field(D);

% restrict to D1-D2
%D = D(1:2);

% sort metrics by combined pval rank
[A, B, rank, idx_A, idx_B, p] = sort_metric_pval_ranks(D,'Trim',false);

% scatter results
scatter_field_pairs(D(1), [idx_A(1:25) idx_B(1:25)]);


