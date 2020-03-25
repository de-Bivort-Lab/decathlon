function [idx,row_subs,col_subs] = upper_triangle_idx(dim,varargin)

mode = 1;
if ~isempty(varargin)
    mode = varargin{1};
end


if dim > 1
    L=1:dim;
    if mode
        subs = arrayfun(@(x) [L(L<x)' repmat(x,sum(L<x),1)],L,'UniformOutput',false);
    else
        subs = arrayfun(@(x) [repmat(x,sum(L<x),1) L(L<x)'],L,'UniformOutput',false);
    end   
    subs = cat(1,subs{:});
    idx = sub2ind([L(end) L(end)],subs(:,1),subs(:,2));
    row_subs = subs(:,1);
    col_subs = subs(:,2);
else
    idx = 1;
    row_subs = 1;
    col_subs = 1;
end