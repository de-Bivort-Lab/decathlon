function D_paired = pair_decathlon_structs(D,varargin)
% Permutes the field order of decathlon data structs in struct array D to match
% and (optionally) collapses metrics in each struct into distilled form via PCA.

% parse processing options and specify default options if none are specified
opts = parse_processing_options(varargin{:});
collapse_mode = opts{find(strcmp(opts,'CollapseMode'))+1};

% match fields between decathlon structs
for i=1:numel(D)
    for j=1:numel(D)-1
        other_idx = find(~ismember(1:numel(D),i));
        [D(i),D(other_idx(j))] = pairFields(D(i),D(other_idx(j)),opts{:});
    end
end

% collapse metrics
if strcmpi(collapse_mode,'PCA')
    D = combine_PCA_split(D,opts{:});
else
    for i=1:numel(D)
        D(i) = collapseMetrics(D(i),opts{:});
    end
end
D_paired = D;

