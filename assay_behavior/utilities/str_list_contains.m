function list_contains = str_list_contains(list,patterns)

% make case insensitive
list = lower(list);
patterns = lower(patterns);
if ~iscell(patterns)
    patterns = {patterns};
end


list_contains = false(numel(list),numel(patterns));
for i=1:numel(patterns)
    list_contains(:,i) = cellfun(@(str) any(strfind(str,patterns{i})), list);
end