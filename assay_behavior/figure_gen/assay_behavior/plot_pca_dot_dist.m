function dot_prods = plot_pca_dot_dist(pcs_a,pcs_b)

% iterate over PCs in A
dot_prods = NaN(size(pcs_a,2));
for i=1:size(pcs_a,2)
    % dot PCs together
    dp = abs(cellfun(@(b) dot(pcs_a(:,i),b), num2cell(pcs_b,1)));
    dot_prods(i,:) = sort(dp,'descend');
end

imagesc(dot_prods);
