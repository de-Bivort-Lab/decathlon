function trans_prob = find_mode_transitions(start_idx,unique_modes,iter,max_iter)
% calculate the transition probabilities from mode A to mode B
% 
% output:   trans_prob - element i,j is the prob. trans from i to j

% print the iteration
fprintf('calculating transition probabilities for fly %i of %i\n',iter,max_iter)

% get sequence of individual behavioral modes
mode_seq = cellfun(@(si,i) repmat(i,numel(si),1), ...
    start_idx, num2cell(1:numel(start_idx))', 'UniformOutput', false);
mode_seq = cat(1,mode_seq{:});
start_idx = cat(1,start_idx{:});

% sort by the order of mode bout starting indices
[~,perm] = sort(start_idx);
mode_seq = mode_seq(perm);

% iterate over each behavioral mode and count the transitions to other modes
max_idx = numel(mode_seq);
mode_start_idx = arrayfun(@(mode_idx) find(mode_seq==mode_idx),...
    unique_modes,'UniformOutput', false);
next_mode_idx = cellfun(@(msi) msi(msi<max_idx)+1,...
    mode_start_idx,'UniformOutput',false);
next_mode_cts = cellfun(@(nmi) histc(mode_seq(nmi)',unique_modes),...
    next_mode_idx,'UniformOutput',false);
trans_prob = cat(1,next_mode_cts{:})';
trans_prob = trans_prob./sum(trans_prob);

