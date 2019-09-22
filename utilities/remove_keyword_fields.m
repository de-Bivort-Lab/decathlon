function D = remove_keyword_fields(D, keywords, varargin)
% remove fields that contain any keywords

kw_grp = '';
for i=1:numel(varargin)
    arg = varargin{i};
    if ischar(arg)
        switch arg
            case 'KeywordGroup'
                i = i+1;
                kw_grp = varargin{i};
        end
    end
end

switch kw_grp
    case 'redundant'
        keywords = {'nTrials';'preodor_occupancy';'nBouts';'Y-Maze circling';...
            'bout_length';'Culling'};
end

for i=1:numel(D)
    f = D(i).fields;
    remove = any(str_list_contains(f,keywords),2);
    fn = fieldnames(D(i));
    for j=1:numel(fn)
        tmp_field = D(i).(fn{j});
        if size(tmp_field,1) == numel(f)
            tmp_field = tmp_field(~remove,:);
        end
        if size(tmp_field,2) == numel(f)
            tmp_field = tmp_field(:,~remove);
        end
        D(i).(fn{j}) = tmp_field;
    end
end
            