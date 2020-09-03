function D = load_decathlon_structs(fdir,var_name)
% loads decathlon data struct matching var_name from decathlon_final_data.mat
% in the parent directory (fdir) 

% define file path
fpaths = recursiveSearch(fdir);
[~,fnames,fext] = cellfun(@fileparts,fpaths,'UniformOutput',false);
fnames = cellfun(@(path,ext) [path ext], fnames, fext, 'UniformOutput', false);
fidx = find(strcmp(fnames,'decathlon_behavior_data.mat'));
if ~fidx
    error('Could not locate decathlon_behavior_data.mat in input directory.');
end

D = load(fpaths{fidx},var_name);
fn = fieldnames(D);
D = D.(fn{1});

