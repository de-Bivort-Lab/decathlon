function D = apply_cluster_permutation(perm, D)

for i=1:numel(D)
    D(i).data = D(i).data(:,perm);
    D(i).fields = D(i).fields(perm);
end