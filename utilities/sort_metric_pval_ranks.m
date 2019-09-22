function [field_A, field_B, ranks, idx_A, idx_B, p_vals] = ...
    sort_metric_pval_ranks(D, varargin)
% sort metric pairs by their combined pval ranks, from lowest to highest

% parse inputs
collapse_mode = 'average';
collapse_fields = 'none';
do_trim = true;
for i=1:numel(varargin)
    arg = varargin{i};
    if ischar(arg)
        switch arg
            case 'CollapseFields'
                i = i+1;
                collapse_fields = varargin{i};
            case 'CollapseMode'
                i = i+1;
                collapse_mode = varargin{i};
            case 'Trim'
                i = i+1;
                do_trim = varargin{i};
        end
    end
end

% initialize linear indices of unique field pairs
pair_idx = upper_triangle_idx(numel(D(1).fields));

% initialize p-value and rank placeholders
p_vals = cell(numel(D),1);
p_rank = cell(numel(D),1);

for i=1:numel(D)
    % calculate covariance matrix
    [r,p] = corr(D(i).data,'Type','spearman','rows','pairwise');

    % replace NaNs
    r(isnan(r))=0;
    p(isnan(p))=1;
    p(r==1) = 1;
    
    % sort the p_values
    pair_p = p(pair_idx);
    [sorted_p, ~] = sort(pair_p);
    p_vals{i} = p;
    ranks = arrayfun(@(pp) find(sorted_p==pp,1), pair_p);
    p_rank{i} = ranks;
end

% sort all metric pairs by the lowest combined rank
p_rank = cat(2,p_rank{:});
t_rank = sum(p_rank,2);
[total_rank, permutation] = sort(t_rank);
ranks = p_rank(permutation,:);

% convert to metric pairs
ranked_idx = pair_idx(permutation);
[idx_A, idx_B] = ind2sub(size(p), ranked_idx);
field_A = D(1).fields(idx_A);
field_B = D(1).fields(idx_B);