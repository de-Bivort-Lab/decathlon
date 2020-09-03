

a = cat(1,D_enrichment_results.cat_genes);
a = cat(1,a{:});
a = unique(a);


%%

% read in raw text
p = ['C:\Users\winsl0w\Documents\decathlon\decathlon_paper_code\'...
    'decathlon\gene_expression\kegg_all_gene_id_to_uniprot.txt'];

fid = fopen(p);
raw = textscan(fid,'%s','delimiter','\t');
raw = raw{1};
fclose(fid);

% parse pathways and genes
k = raw(logical(mod(1:numel(raw),2)));
k = regexp(k,'(?<=:).*','match');
k = cellfun(@(s) s{1}, k, 'UniformOutput', false);
u = raw(~mod(1:numel(raw),2));
u = regexp(u,'(?<=:).*','match');
u = cellfun(@(s) s{1}, u, 'UniformOutput', false);

% filter out duplicate IDs
[k,pk] = unique(k,'stable');
u = u(pk);

%%

p_out = ['C:\Users\winsl0w\Documents\decathlon\decathlon_paper_code\'...
    'decathlon\gene_expression\all_uniprot.txt'];
fid = fopen(p_out,'W+');
cellfun(@(s) fprintf(fid,'%s\n',s), u);
fclose(fid);

%%

% read in raw text
p = ['C:\Users\winsl0w\Documents\decathlon\decathlon_paper_code\'...
    'decathlon\gene_expression\kegg_all_gene_id_to_ncbi_id.txt'];

fid = fopen(p);
raw = textscan(fid,'%s','delimiter','\t');
raw = raw{1};
fclose(fid);

% parse pathways and genes
k = raw(logical(mod(1:numel(raw),2)));
k = regexp(k,'(?<=:).*','match');
k = cellfun(@(s) s{1}, k, 'UniformOutput', false);
u = raw(~mod(1:numel(raw),2));
u = regexp(u,'(?<=:).*','match');
u = cellfun(@(s) s{1}, u, 'UniformOutput', false);

% filter out duplicate IDs
[k,pk] = unique(k,'stable');
u = u(pk);

%%

p_out = ['C:\Users\winsl0w\Documents\decathlon\decathlon_paper_code\'...
    'decathlon\gene_expression\all_ncbi_id.txt'];
fid = fopen(p_out,'W+');
cellfun(@(s) fprintf(fid,'%s\n',s), u);
fclose(fid);

%% parse kegg-to-ncbi list for unique FBgn

p_in = ['C:\Users\winsl0w\Documents\decathlon\decathlon_paper_code\'...
    'decathlon\gene_expression\all_ncbi_gid_2_fbgn_lookup.txt'];
fid = fopen(p_in);
raw = textscan(fid,'%s','delimiter','\t');
raw = raw{1};
fclose(fid);

% parse pathways and genes
ncbi = raw(mod(1:numel(raw),3)==1);
fbgn = raw(mod(1:numel(raw),3)==2);
sym = raw(mod(1:numel(raw),3)==0);
all_ids = [ncbi(2:end) fbgn(2:end) sym(2:end)];

% remove clone gene ids
is_fbgn = ~cellfun(@isempty,regexp(all_ids(:,2),'FBgn','match'));
all_ids = all_ids(is_fbgn,:);

% find duplicate ncbi-ds
[u_ncbi,~,p_ncbi] = unique(all_ids(:,1),'stable');
hc = histc(p_ncbi,1:numel(u_ncbi));
dup_idx = find(hc>1);
dup_pair_idx = cell(size(dup_idx));
for i=1:numel(dup_idx)
   dup_pair_idx{i} = find(strcmp(all_ids(:,1),u_ncbi(dup_idx(i))));
end

p_out = ['C:\Users\winsl0w\Documents\decathlon\decathlon_paper_code\'...
    'decathlon\gene_expression\all_kegg_fbgn_query.txt'];
fid = fopen(p_out,'W+');
cellfun(@(s) fprintf(fid,'%s\n',s),all_ids(:,2));
fclose(fid);



%%

f1 = lookup_table.fbgn;
f2 = geneid_lookup_tbl.fbgn;
f = [f2;f1];

k1 = lookup_table.kegg;
k2 = geneid_lookup_tbl.kegg;
n_kegg = numel(unique([k1;k2]));
k = [k2;k1];

% remove invalid kegg ids
valid_kegg = ~cellfun(@isempty,regexp(k,'Dmel_','match'));
k = k(valid_kegg);
f = f(valid_kegg);
[uk,pk] = unique(k);
k = uk;
f = f(pk);

% remove invalid FBgn
valid_fbgn = ~cellfun(@isempty,regexp(f,'FBgn','match'));
k = k(valid_fbgn);
f = f(valid_fbgn);
[uf,pf] = unique(f);
f = uf;
k = k(pf);

k(strcmp(f,'FBgn0037417')) = [];
f(strcmp(f,'FBgn0037417')) = [];
k(strcmp(f,'FBgn0036125')) = [];
f(strcmp(f,'FBgn0036125')) = [];

p_out = ['C:\Users\winsl0w\Documents\decathlon\decathlon_paper_code\'...
    'decathlon\gene_expression\all_kegg_fbgn_query.txt'];
fid = fopen(p_out,'W+');
cellfun(@(s) fprintf(fid,'%s\n',s),f);
fclose(fid);

%% after submitting exported query list of FBgn, parse the queried field file

p_in = ['C:\Users\winsl0w\Documents\decathlon\decathlon_paper_code\'...
    'decathlon\gene_expression\all_kegg_fbgn_field_file.txt'];
fid = fopen(p_in);
raw = textscan(fid,'%s','delimiter','\t');
raw = raw{1};
fclose(fid);

% parse pathways and genes
ncols = 4;
validated_fbgn = raw(mod(1:numel(raw),ncols)==2);
gene_name = raw(mod(1:numel(raw),ncols)==3);
gene_symbol = raw(mod(1:numel(raw),ncols)==0);
field_file = [validated_fbgn(2:end) gene_symbol(2:end) gene_name(2:end)];

% match kegg-ncbi lookup table to field file
pk = cellfun(@(s) find(strcmp(f,s)), field_file(:,1));
field_file = [field_file k(pk)];
field_file = field_file(:,[1 4 2 3]);
lookup_table = table(field_file(:,1),field_file(:,2),...
    field_file(:,3),field_file(:,4),'VariableNames',{'fbgn';'kegg';'symbol';'name'});

