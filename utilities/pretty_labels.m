function fields = pretty_labels(fields)

for i=1:numel(fields)
    tmp = fields{i};
    tmp(tmp=='_') = ' ';
    fields{i} = tmp;
end