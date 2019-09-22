function n = get_pairwise_sampling(data)

nFields = size(data,2);
n = NaN(nFields);
for i=1:nFields
    for j=1:nFields       
        if i~=j
            n(i,j) = sum(~isnan(data(:,i)) & ~isnan(data(:,j)));
        end
    end
end