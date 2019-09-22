function D = standardize_by_field(D)

for i=1:numel(D)
   D(i).data = nanzscore(D(i).data); 
end