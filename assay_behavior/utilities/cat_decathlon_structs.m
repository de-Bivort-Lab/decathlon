function D_cat = cat_decathlon_structs(D,varargin)

% parse inputs
opts = parse_processing_options(varargin{:});

% collapse metrics
for i=1:numel(D)
    D(i) = collapseMetrics(D(i),opts{:});
end

% pair the matrices
D_p = pair_decathlon_structs(D,opts{:});
paired_fields = standardize_fieldnames(D_p(1).fields);
[assay, metric] = parse_fieldnames(paired_fields);
paired_fields = cellfun(@(a,m) sprintf('%s %s',a,m), assay, metric,...
    'UniformOutput', false);
unpaired_fields = cell(numel(D),1);
for i = 1:numel(D)
    f = standardize_fieldnames(D(i).fields);
    [assay, metric] = parse_fieldnames(f);
    f = cellfun(@(a,m) sprintf('%s %s',a,m), assay, metric,...
        'UniformOutput', false);
    unpaired_fields{i} = find(~any(str_list_contains(f,paired_fields),2));
end

% concatenate the matrices
n = arrayfun(@(d) size(d.data,1), D);
if size(n,1) > size(n,2)
    n = n';
end
n = [0 cumsum(n)];
D_cat.data = cat(1,D_p.data);
if isfield(D_p(1),'imputed')
    D_cat.imputed = cat(1,D_p.imputed);
end
D_cat.fields = D_p(1).fields;
f = fieldnames(D_p(1).meta);
for i=1:numel(f)
    new_meta = arrayfun(@(d) d.meta.(f{i}), D_p, 'UniformOutput', false);
    D_cat.meta.(f{i}) = cat(1,new_meta{:});
end
for i=1:numel(D)
    new_data = D(i).data(:,unpaired_fields{i});
    unpaired_data = NaN(size(D_cat.data,1),size(new_data,2));
    unpaired_data(n(i)+1:n(i+1),:) = new_data;
    D_cat.data =  [D_cat.data unpaired_data];
    D_cat.fields = [D_cat.fields; D(i).fields(unpaired_fields{i})];
    if isfield(D(i),'imputed')
        new_imputed = D(i).imputed(:,unpaired_fields{i});
        unpaired_imputed = true(size(D_cat.data,1),size(new_imputed,2));
        unpaired_imputed(n(i)+1:n(i+1),:) = new_imputed;
        D_cat.imputed =  [D_cat.imputed unpaired_imputed];
    end
    
    f = fieldnames(D_cat.meta);
    for j=1:numel(f)
        new_meta = D(i).meta.(f{j})(:,unpaired_fields{i});
        if iscell(new_meta)
            unpaired_meta = cell(size(D_cat.data,1),size(new_meta,2));
            unpaired_meta(n(i)+1:n(i+1),:) = new_meta;
        else
            unpaired_meta = NaN(size(D_cat.data,1),size(new_meta,2));
            unpaired_meta(n(i)+1:n(i+1),:) = new_meta;
        end
        D_cat.meta.(f{j}) =  [D_cat.meta.(f{j}) unpaired_meta];
    end
end

D_cat.n = get_pairwise_sampling(D_cat.data);
