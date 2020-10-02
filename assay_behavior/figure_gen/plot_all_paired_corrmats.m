% define a save path
impute_mode = 'als_avg';
save_dir = uigetdir('D:/decathlon_data_and_analysis/decathlon_analysis/figures/',...
    'Select path to save figure outputs');
name = sprintf('D12_D3_batch_standardized_bout_reprocessed_%s_filled',impute_mode);
ext = '.fig';

D = impute_decathlon_structs(D,'ImputeMode',impute_mode);
D = standardize_by_field(D);

[fh1, fh2] = plot_paired_corrmats(D,'CollapseFields','none');
saveas(fh1,cat(2,save_dir,'\',name,'corrmat_raw',ext));

[fh1, fh2] = plot_paired_corrmats(D,'CollapseFields','circadian');
saveas(fh1,cat(2,save_dir,'\',name,'corrmat_circadian-collapsed',ext));

[fh1, fh2] = plot_paired_corrmats(D,'CollapseFields','all','Cluster',true);
saveas(fh1,cat(2,save_dir,'\',name,'corrmat_apriori-collapsed',ext));

[fh1, fh2] = plot_paired_corrmats(D,'CollapseFields','all','Cluster',true,...
    'CollapseMode','PCA','PCs',2,'Standardize',false,'ImputeMode','none');
saveas(fh1,cat(2,save_dir,'\',name,'corrmat_apriori-PCA',ext));