function [apriori_grouped, apriori_group_names, grp_idx] = group_apriori_pcs(D)

f = D.fields;
grp = cellfun(@(s) s(1:find(s=='(',1)-2), f, 'UniformOutput', false);
[apriori_group_names,~,c] = unique(grp,'stable');

ngrps = numel(apriori_group_names);
apriori_grouped = cell(ngrps,1);
grp_idx = cell(ngrps,1);
for i=1:ngrps
    grp_idx{i} = find(c==i); 
    apriori_grouped{i} = D.data(:,c==i);
end