function D = standardize_data_by_group(D)

for i=1:numel(D)
    data = D(i).data;
    for j = 1:size(data,2)
       timeofday = D(i).meta.timeofday(:,j);
       grps = unique(timeofday);
       metric = data(:,j);
       for k=1:numel(grps)
          mask = timeofday == grps(k);
          grp_data = metric(mask);
          metric(mask) = nanzscore(grp_data);
       end
       data(:,j) = metric;
    end
    D(i).data = data;
end