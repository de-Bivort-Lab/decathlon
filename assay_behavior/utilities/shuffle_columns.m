function data = shuffle_columns(data)

[nr nc] = size(data);
for i=1:nc
   data(:,i) = data(randperm(nr),i);
end