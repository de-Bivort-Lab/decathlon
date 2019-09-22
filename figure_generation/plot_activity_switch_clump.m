
% initialize struct array for paired data
D_p = D;

% match fields between decathlon structs
for i=1:numel(D_p)
    for j=1:numel(D_p)-1
        other_idx = find(~ismember(1:numel(D_p),i));
        [D_p(i),D_p(other_idx(j))] = pairFields(D_p(i),D_p(other_idx(j)),'Trim',true);
    end
end

% re-arrange in clustered order
p = get_cluster_permutation(D_p);
D_p = apply_cluster_permutation(p, D_p);


figure;
for i=1:numel(D_p)
    ah1 = subplot(2,numel(D_p),i);
    ah2 = subplot(2,numel(D_p),numel(D_p)+i);
   corrMat_byGroup({{'switchiness';'clumpiness';'speed';'nTrials'}},D_p(i),...
       'PlotCorrOptions',{'Cluster',false,'Patch',false,'Parent',[ah1 ah2],'Signed',false}); 
end