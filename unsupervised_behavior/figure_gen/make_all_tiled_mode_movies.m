function make_all_tiled_mode_movies(movieVec,aligned_movie_dir,L)


vidFiles = recursiveSearch(aligned_movie_dir,'ext','.avi');

if length(vidFiles) > L

    vidFiles = reshape(vidFiles',4,length(vidFiles)/4)';
    [~,fLabels,~] = cellfun(@fileparts,vidFiles(:,1),'UniformOutput',false);
    str_idx = cellfun(@(fL) find(fL=='_'),fLabels,'UniformOutput',false);
    strains = cellfun(@(fL,sidx) fL(sidx(2)+1:sidx(3)-1),...
        fLabels,str_idx,'UniformOutput',false);
    id_nums = cellfun(@(fL,sidx) str2double(fL(sidx(3)+1:sidx(4)-1)),fLabels,str_idx);
    vidIDs = cellfun(@(s,id) sprintf('%s_%03.0f',s,id),strains,num2cell(id_nums),...
        'UniformOutput', false);
    [~,p] = sort(vidIDs);
    vidFiles = num2cell(vidFiles,2);
    vidFiles = vidFiles(p);
end

vidFiles = cat(1,vidFiles{:});
nFrames = get_movie_framenumbers(vidFiles);

% prompt user for save path
saveDir = uigetdir('Select directory to save output movie files');
nTiles = [8 8];
for i=1:size(movieVec,1)
    savePath = [saveDir '\tiledSample_mode' num2str(i) '.mp4'];
    modeMovie(movieVec{i},vidFiles,nTiles,nFrames,100,savePath,i);
end