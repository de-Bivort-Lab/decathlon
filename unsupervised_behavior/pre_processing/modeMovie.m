function modeMovie(movieVecs,vidPaths,nTiles,nFrames,dur,savePath,modeNum)

% restrict the movie vector to the first nVids unique flies
nVids = prod(nTiles);

% sort by strain
strain_id = movieVecs(:,2) < 168;
vec_idx = NaN(nVids/2,2);

for strain=0:1
    i=1; j=1;
    strain_idx = strain_id==strain;
    tmp_vec = movieVecs(strain_idx,:);
    unique_ids = NaN(nVids,1);
    while any(isnan(vec_idx(:,strain+1))) && i < size(tmp_vec,1)
        k = tmp_vec(i,2);
        if all(~ismember(k,unique_ids))
           vec_idx(j,strain+1) = i;
           unique_ids(j) = k;
           j = j+1;
        end
        i=i+1;
    end
    strain_idx = find(strain_idx);
    tmp_idx = vec_idx(~isnan(vec_idx(:,strain+1)),strain+1);
    vec_idx(1:numel(tmp_idx),strain+1) = strain_idx(tmp_idx);
end
vec_idx(isnan(vec_idx)) = [];
movieVecs = movieVecs(vec_idx(:),:);

nBout = size(movieVecs,1);
dur = min(movieVecs(:,3));

q = VideoReader(vidPaths{1,1});
xPix = q.Width;
yPix = q.Height;
delete(q);

frameRange = cumsum(nFrames,2);
newVid = VideoWriter(savePath,'MPEG-4');
newVid.FrameRate = 20;
open(newVid);

fprintf(1,'\t Initializing mode #%2i video\n',modeNum);

% initialize videoReader objects
vid_objs = VideoReader.empty(nVids,0);
vid_idx = NaN(nVids,2);
for i = 1:nTiles(1)
    for j = 1:nTiles(2)
        listIdx = (i-1)*nTiles(1)+j;
        
        if listIdx <= nBout
            frame_idx = movieVecs(listIdx,1)+1;
            vid_idx(listIdx,1) = movieVecs(listIdx,2);
            vid_idx(listIdx,2) = find(frame_idx <= frameRange(vid_idx(listIdx,1),:),1,'first');
            fprintf('initializing videoReader %i of %i\n',listIdx,nVids);
            vid_objs(listIdx) = VideoReader(vidPaths{vid_idx(listIdx,1),vid_idx(listIdx,2)});
        end
    end
end

%%
for k = 1:dur
    
    fprintf(1,'\t Writing frame #%2i of %2i\n',k,dur);
    
    fIm = uint8(zeros(xPix*nTiles(1),yPix*nTiles(2),3));
    
    for i = 1:nTiles(1)
        for j = 1:nTiles(2)
            
            listIdx = (i-1)*nTiles(1)+j;
            if listIdx <= nBout

                % get current frame number
                fn = movieVecs(listIdx,1)+k;
                
                % update video object if necessary
                if fn > frameRange(vid_idx(listIdx,1),vid_idx(listIdx,2))
                    delete(vid_objs(listIdx));
                    vid_idx(listIdx,2) = vid_idx(listIdx,2)+1;
                    vid_objs(listIdx) = VideoReader(vidPaths{vid_idx(listIdx,1),vid_idx(listIdx,2)});
                end
                
                % adjust frame number relative to subvideo
                if vid_idx(listIdx,2) > 1
                    fn = fn - frameRange(vid_idx(listIdx,1),vid_idx(listIdx,2)-1);
                end

                % read image for current tile
                subIm = read(vid_objs(listIdx),fn);
                subIm=subIm(:,:,1);
                subIm = repmat(subIm,1,1,3);
                subIm(10:20,10:20,:) = 0;
                if listIdx > nVids/2
                    subIm(10:20,10:20,3) = 255;
                else
                    subIm(10:20,10:20,1) = 255;
                end
                
                fIm((j-1)*yPix+1:j*yPix,(i-1)*xPix+1:i*xPix,:) = subIm;
            end   
        end
    end
    
    writeVideo(newVid,fIm);
    
end

close(newVid);