function pairs = unique_idx_pairs(n,mode)

pairs = NaN(sum(1:n),2);
pct = 0;
for i=1:n
    for j=1:n
        if j >= i+mode
            pct = pct+1;
            pairs(pct,:) = [i j];
        end
    end
end

pairs(isnan(pairs(:,1)),:) = [];