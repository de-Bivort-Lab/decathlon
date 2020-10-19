function data_struct = load_decathlon_structs(fdir,var_name)
% loads decathlon data struct matching var_name from decathlon_final_data.mat
% in the parent directory (fdir) 

% define file path
fpaths = recursiveSearch(fdir);
data_struct = [];

fnames = {'decathlon_assay_data';'decathlon_enrichment_results';...
    'decathlon_rnaseq_data';'decathlon_thermo_gal4_data';...
    'decathlon_unsupervised_behavior_data';'decathlon_unsupervised_embedding_timeseries'};
[~,all_fnames,~] = cellfun(@(s) fileparts(s), fpaths, 'UniformOutput', false);
var_names = {{'D_als_filled';'D_als_filled_batch_merged';'D_raw_unfilled';'D_zscored_unfilled'};...
                {'D_enrichment_results';'D_rnaseq_models'};...
                {'D_seq'};
                {'D_thermo'};...
                {'D_us'};...
                {'embeddings'}};

has_match = cellfun(@(s) any(strcmpi(s,var_name)), var_names);
if any(has_match)
   fidx = find(strcmpi(all_fnames,fnames{has_match}),1);
   D = load(fpaths{fidx});
   fn = fieldnames(D);
   data_struct = D.(fn{strcmpi(fn,var_name)});
end

if isempty(data_struct)
   error('No data struct with that name found in the input directory.'); 
end



