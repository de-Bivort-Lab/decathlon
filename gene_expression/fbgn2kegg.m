function kegg_gene_list = fbgn2kegg(fbgn_gene_list,lookup_table)
% convert flybase gene no. to kegg gene ID

% find indices in lookup table
lookup_idx = cellfun(@(g) find(strcmp(lookup_table.fbgn,g)),...
    fbgn_gene_list, 'UniformOutput', false);
lookup_idx = cat(1,lookup_idx{:});
kegg_gene_list = lookup_table.kegg(lookup_idx);

% format kegg ids
kegg_gene_list = cellfun(@(s) strsplit(s,'; '),kegg_gene_list, 'UniformOutput', false);
kegg_gene_list = cat(2,kegg_gene_list{:})';
kegg_gene_list(strcmp(kegg_gene_list,'-')) = [];
if iscell(kegg_gene_list)
    is_melanogaster = cellfun(@(s) strcmp(s(1:3),'dme'), kegg_gene_list);
    kegg_gene_list(~is_melanogaster) = [];
    kegg_gene_list = cellfun(@(s) s(5:end), kegg_gene_list, 'UniformOutput', false);
else
    kegg_gene_list = {};
end
