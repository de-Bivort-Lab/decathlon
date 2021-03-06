function [grouped_fields, group_idx] = groupFields(fields, keywords)


if ~iscell(keywords) && ischar(keywords)
   keywords = {keywords};
end

if any(~cellfun(@(f) ischar(f), fields))
   error('All keywords must be a string'); 
end

include = false(numel(keywords),numel(fields));
for i=1:numel(keywords)
    include(i,:) = any(str_list_contains(fields,keywords{i}),2);
end

grouped_fields = fields(any(include,1));
group_idx = find(any(include,1));