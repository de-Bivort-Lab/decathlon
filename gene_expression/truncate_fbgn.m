function int = truncate_fbgn(fbgn)
    int = regexp(fbgn,'(?<=(FBgn0)).*','match');
    int = cat(1,int{:});
    if iscell(int)
        int = cellfun(@str2double,int);
    else
        int = str2double(int);
    end
    int = uint32(int);
end