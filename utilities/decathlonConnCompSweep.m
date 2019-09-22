function [out,comps,all_bins]=decathlonConnCompSweep(data,reps)

numDims=size(data,2);

data(all(isnan(data),2),:) = [];
cov_mat = abs(cov(zscore(data)));
cov_vals = cov_mat(upper_triangle_idx(size(cov_mat,1)));
lower_bound = prctile(cov_vals,1);
upper_bound = max(cov_mat(:));
%upper_bound = max(cov_vals);
%cov_mat(logical(diag(ones(size(cov_mat,1),1)))) = max(cov_vals(:)).*1.01;
vals=linspace(0,1,reps);

comps=[];
all_bins = cell(reps,1);
for i=1:length(vals)
    G=graph(cov_mat>vals(i));
    bins=conncomp(G);
    comps=[comps;length(unique(bins))];
    all_bins{i} = bins;
end

out=hist(comps,1:numDims);
%out(1) = 0;
out = out./sum(out);