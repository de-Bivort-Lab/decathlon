function movieVector = getModeMovieVector(frameIdx,durations)

% generate an N x 2 vector where the first column is the first frame number,
% the second is fly number, and the third is duration
nVids = numel(frameIdx);
nSamplesPerVid = cellfun(@numel,frameIdx);
idList = arrayfun(@makeIdList,nSamplesPerVid,1:nVids,'UniformOutput',false);
idList = cat(1,idList{:});
frameIdx = cat(1,frameIdx{:});
durations = cat(1,durations{:});
[~,p] = sort(durations);
p = fliplr(p')';

movieVector = [frameIdx(p) idList(p) durations(p)];


function out = makeIdList(n,idx)

out = ones(n,1).*idx;