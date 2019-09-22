function status = isfield_all(s,field)

fields = regexp(field,'\.','split');
try
    status = true;
    for i=1:numel(fields)
        status = (isfield(s,fields{i}) || isprop(s,fields{i})) && status;
        s = s.(fields{i});
    end
catch
    status = false;
end