function turns = parse_v2_ymaze_data(flyTracks)

% calculate speed
cen = flyTracks.centroid;
cen_previous = [NaN(1,2,size(cen,3)); flyTracks.centroid(1:size(cen,1)-1,:,:)];
cen_current = flyTracks.centroid;
cen_next = [flyTracks.centroid(2:size(cen,1),:,:); NaN(1,2,size(cen,3))];
cen = cat(4,cen_previous,cen_current,cen_next);
cen = nanmean(cen,4);
speed = calculate_speed(cen,flyTracks.tStamps);
[idx,lag_thresh,spd_thresh,fh] = blockActivity_v2(speed);

% record bout metrics
idx(cellfun(@isempty,idx)) = {zeros(0,2)};
bout_temp = cellfun(@(i) ceil(flyTracks.tStamps(i(:,1))./flyTracks.tStamps(end)),...
    idx,'UniformOutput', false);
nBouts = NaN(size(speed,2),1);
bout_length = NaN(size(speed,2),1);
bout_clumpiness = NaN(size(speed,2),1);
bout_idx = cell(size(speed,2),1);

t = flyTracks.tStamps;
speed_blocks = cell(size(nBouts));
mean_speed = NaN(size(nBouts));
temp_idx = ceil(t./3600);
for i=1:1
    bout_idx(:,i) = cellfun(@(ii,bt) ii(bt==i,:), idx, bout_temp, 'UniformOutput', false);
    nBouts(:,i) = cellfun(@(i) size(i,1), bout_idx(:,i));
    bout_length(:,i) = cellfun(@(ii) nanmean(diff(ii,1,2)),bout_idx(:,i));
    
    % calculate bout length clumpiness
    dur = t(end)/1;
    iti = cellfun(@(ii) diff([t(ii(1:end-1,2)) t(ii(2:end,1))],1,2), bout_idx(:,i), 'UniformOutput', false);
    bout_clumpiness(:,i) = cellfun(@(i,n) std(i)/(dur/n), iti, num2cell(nBouts(:,i)));
    
    % speed measures
    speed_blocks(:,i) = num2cell(speed(temp_idx==i,:),1);
    mean_speed(:,i) = nanmean(speed(temp_idx==i,:));
end

% get turn indices
turns.idx = ~isnan(flyTracks.rightTurns);
turns.idx = num2cell(turns.idx,1);
turns.idx = cellfun(@find,turns.idx,'UniformOutput',false);

% convert sequence of maze arms to left/right turn sequence
arm_sequence = cellfun(@(aseq,ti) aseq(ti), ...
    num2cell(flyTracks.rightTurns,1), turns.idx, 'UniformOutput', false);
turn_sequence = cellfun(@score_turns,arm_sequence,num2cell(flyTracks.mazeOri)',...
    'UniformOutput', false);

% remove first turn
no_turns = cellfun(@numel,turns.idx)<1;
trim_turns = cellfun(@(ti) ti(2:end), turns.idx(~no_turns), 'UniformOutput', false);
turns.idx(no_turns) = {[]};
turns.idx(~no_turns) = trim_turns;
turns.n = cellfun(@numel,turns.idx);
clear trim_turns no_turns

% Calculate turn metrics
turns.t = cellfun(@(ti) flyTracks.tStamps(ti), turns.idx, 'UniformOutput', false);
turns.right_bias = cellfun(@sum,turn_sequence)./turns.n;
turns.switchiness = cellfun(@(s,r,nt) sum((s(1:end-1)+s(2:end))==1)/(2*r*(1-r)*nt),...
    turn_sequence,num2cell(turns.right_bias),num2cell(turns.n));
turns.clumpiness = cellfun(@(t,nt) std(diff([0;t]))/(flyTracks.tStamps(end)/nt),...
    turns.t, num2cell(turns.n));

% calculate LED ymaze metrics if applicable
if strcmpi(flyTracks.exp,'LED Y-maze')
    light_sequence = cellfun(@(lseq) lseq(~isnan(lseq)), ...
        num2cell(flyTracks.lSeq,1), 'UniformOutput', false);
    turns.light_bias = cellfun(@sum,light_sequence)./(turns.n+1);
    turns.light_switchiness = cellfun(@(s,r,nt) sum((s(1:end-1)+s(2:end))==1)/(2*r*(1-r)*nt),...
        turn_sequence,num2cell(turns.right_bias),num2cell(turns.n+1));
end

% store metrics in struct
turns.bout_clumpiness = bout_clumpiness;
turns.nBouts = nBouts;
turns.bout_length = bout_length;
turns.speed = mean_speed;

fn = fieldnames(turns);
for i=1:numel(fn)
    tmp = turns.(fn{i});
    if any(size(tmp)==flyTracks.nFlies) && size(tmp,1)~=flyTracks.nFlies
        tmp = tmp';
    end
    turns.(fn{i}) = tmp;
end



function speed = calculate_speed(cen,t)

dist = [NaN(1,size(cen,3));squeeze(sqrt(sum(diff(cen,1,1).^2,2)))];
speed = dist./[NaN;diff(t)];


function turn_seq = score_turns(arm_seq,ori)

tSeq=diff(arm_seq);  
if ori
    turn_seq = tSeq==1 | tSeq==-2;
else
    turn_seq = tSeq==-1 | tSeq==2;
end

