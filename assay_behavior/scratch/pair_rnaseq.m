function D = pair_rnaseq(D)

pairs = unique_idx_pairs(numel(D),1);
for i=1:size(pairs,1)
    % restrict geneIDs and data to common genes
    is_common_a = ismember(D(pairs(i,1)).geneID,D(pairs(i,2)).geneID);
    is_common_b = ismember(D(pairs(i,2)).geneID,D(pairs(i,1)).geneID);
    D(pairs(i,1)).geneID = D(pairs(i,1)).geneID(is_common_a);
    D(pairs(i,1)).data = D(pairs(i,1)).data(:,is_common_a);
    D(pairs(i,2)).geneID = D(pairs(i,2)).geneID(is_common_b);
    D(pairs(i,2)).data = D(pairs(i,2)).data(:,is_common_b);
    
    % sort remaining data
    [D(pairs(i,1)).geneID,p_a] = sort(D(pairs(i,1)).geneID);
    D(pairs(i,1)).data = D(pairs(i,1)).data(:,p_a);
    [D(pairs(i,2)).geneID,p_b] = sort(D(pairs(i,2)).geneID);
    D(pairs(i,1)).data = D(pairs(i,1)).data(:,p_b);
end