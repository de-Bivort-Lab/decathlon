function [output_list,has_match] = convert_gene_ids(input_list,input_type,output_type)
% convert flybase gene no. to kegg gene ID

input_type = lower(input_type);
output_type = lower(output_type);
% locate the lookup table file
[src_dir,~,~] = fileparts(mfilename('fullpath'));
load(cat(2,src_dir,'/lookup_table.mat'),'lookup_table');

valid_types = [lookup_table.Properties.VariableNames 'all'];
if ~any(strcmp(valid_types,input_type)) || ~any(strcmp(valid_types,output_type))
        error('Invalid input or output type');
end

if isempty(input_list)
   input_list = lookup_table.(input_type); 
end

% find indices in lookup table
lookup_idx = cellfun(@(g) find(strcmp(lookup_table.(input_type),g)),...
    input_list, 'UniformOutput', false);
has_match = ~cellfun(@isempty,lookup_idx);
lookup_idx = cat(1,lookup_idx{:});

if strcmp(output_type,'all')
    output_list = lookup_table(lookup_idx,:);
else
    output_list = lookup_table.(output_type)(lookup_idx);
end

