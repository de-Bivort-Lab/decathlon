% get file paths for centroid text files and 
cen_dir = {'E:\decathlon unsupervised videos\decathlon 2-2018 raw video';'F:\decathlon 1-2019 raw'};
cen_paths = recursiveSearch(cen_dir,'ext','.txt','keyword','CentroidData');
vid_paths = recursiveSearch(cen_dir,'ext','.avi');
[~,vid_names] = cellfun(@fileparts,vid_paths,'UniformOutput',false);
[par_dirs,fnames] = cellfun(@fileparts,cen_paths,'UniformOutput',false);
[~,par_dirs] = cellfun(@fileparts,par_dirs,'UniformOutput',false);
%%
% initialize data placeholders
strains = cell(numel(vid_paths),1);
IDs = cell(numel(vid_paths),1);
centroids = cell(numel(vid_paths),1);
ct = 0;

for i=1:numel(cen_paths)
    
    fprintf('parsing centroid file %i of %i\n',i,numel(cen_paths));
    
    % use file name to parse out strain, id, and cam number for centroids
    tmp_vids = vid_names(contains(vid_names,par_dirs{i}))';
    idx = cellfun(@(s) find(s=='_'), tmp_vids, 'UniformOutput', false);
    tmp_strains = cellfun(@(s,j) s(j(2)+1:j(3)-1), tmp_vids, idx, 'UniformOutput', false);
    tmp_IDs = cellfun(@(s,j) str2double(s(j(3)+1:j(4)-1)), tmp_vids, idx, 'UniformOutput', false);
    cam_num = cellfun(@(s) str2double(s(end)), tmp_vids);
   
    % permute by camera number
    [~,cam_perm] = sort(cam_num);
    tmp_strains = tmp_strains(cam_perm);
    tmp_IDs = tmp_IDs(cam_perm);
   
    % assign to struct
    strains(ct+1:ct+numel(tmp_strains)) = tmp_strains;
    IDs(ct+1:ct+numel(tmp_strains)) = tmp_IDs;
    strains(ct+1:ct+numel(tmp_strains)) = tmp_strains;

    % read centroid data from file
    cen = dlmread(cen_paths{i});
    cen = permute(reshape(cen,3,4,numel(cen)/12),[3 1 2]);
    cen = cen(:,2:3,:);
    cen = cen(:,:,1:numel(cam_num));
    centroids(ct+1:ct+numel(tmp_strains)) = num2cell(cen,[1 2]);
    
    % update counter
    ct = ct+numel(cam_perm);
end

% data to struct for output to file
unsupervised_cen.data = centroids;
unsupervised_cen.strain = strains;
unsupervised_cen.IDs = IDs;

% remove empty elements
f = cellfun(@isempty,centroids);
unsupervised_cen.data(f) = [];
unsupervised_cen.strain(f) = [];
unsupervised_cen.IDs(f) = [];
