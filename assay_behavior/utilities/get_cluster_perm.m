function [row_perm,col_perm] = get_cluster_perm(data,method,metric)

Z=linkage(data,method,metric);
f=figure;
[~, ~, row_perm]=dendrogram(Z,0);
close(f);
Z=linkage(data',method,metric);
f=figure;
[~, ~, col_perm]=dendrogram(Z,0);
close(f);