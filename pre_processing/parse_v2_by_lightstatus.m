function [out,fh] = parse_v2_by_lightstatus(flyTracks)

% calculate speed
%f = mod(1:size(flyTracks.centroid,1),3)==0;
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
bout_temp = cellfun(@(i) ceil(flyTracks.tStamps(i(:,1))./3600), idx,'UniformOutput', false);
nBouts = NaN(size(speed,2),4);
bout_length = NaN(size(speed,2),4);
bout_clumpiness = NaN(size(speed,2),4);
bout_idx = cell(size(speed,2),4);

t = flyTracks.tStamps;
speed_blocks = cell(size(nBouts));
mean_speed = NaN(size(nBouts));
temp_idx = ceil(t./3600);
for i=1:4
    bout_idx(:,i) = cellfun(@(ii,bt) ii(bt==i,:), idx, bout_temp, 'UniformOutput', false);
    nBouts(:,i) = cellfun(@(i) size(i,1), bout_idx(:,i));
    bout_length(:,i) = cellfun(@(ii) nanmean(diff(ii,1,2)),bout_idx(:,i));
    
    % calculate bout length clumpiness
    dur = t(end)/4;
    iti = cellfun(@(ii) diff([t(ii(1:end-1,2)) t(ii(2:end,1))],1,2), bout_idx(:,i), 'UniformOutput', false);
    bout_clumpiness(:,i) = cellfun(@(i,n) std(i)/(dur/n), iti, num2cell(nBouts(:,i)));
    
    % speed measures
    speed_blocks(:,i) = num2cell(speed(temp_idx==i,:),1);
    mean_speed(:,i) = nanmean(speed(temp_idx==i,:));
end

% read in stim block parsing file and define times that the stimulus switches
stim_blocks = dlmread('E:\Kyobi data\2015.11.04_LightVDark\scripts\lightSequenceScreen.txt');
stim_t = cumsum([0;stim_blocks(:,1)],1);
stim_trans = [1; diff(stim_blocks(:,2))];
[light_start, dark_start] = deal(stim_trans==1,stim_trans==-1);
[t_light_start,t_dark_start] = deal(stim_t(light_start),stim_t(dark_start));

% define starting stimulus and truncate if necessary to match block number
is_light_first = t_light_start(1) < t_dark_start(1);
if numel(t_light_start) > numel(t_dark_start)
    t_light_start(end) = [];
elseif numel(t_dark_start) > numel(t_light_start)
    t_dark_start(end) = [];
end

light = struct('turn_dat',[],'t',[],'temp_block',[]);
dark = struct('turn_dat',[],'t',[],'temp_block',[]);
if is_light_first
    [light.turn_dat,light.t,light.temp_block] = arrayfun(@(tl,td)...
        block_turn_data(flyTracks.tStamps,flyTracks.rightTurns,tl,td),...
        t_light_start,t_dark_start, 'UniformOutput', false);
    [dark.turn_dat,dark.t,dark.temp_block] = arrayfun(@(tl,td)...
        block_turn_data(flyTracks.tStamps,flyTracks.rightTurns,tl,td),...
        t_dark_start,[t_light_start(2:end); flyTracks.tStamps(end)], 'UniformOutput', false);
else
    [dark.turn_dat,dark.t,dark.temp_block] = arrayfun(@(tl,td)...
        block_turn_data(flyTracks.tStamps,flyTracks.rightTurns,tl,td),...
        t_dark_start,t_light_start, 'UniformOutput', false);
    [light.turn_dat,light.t,light.temp_block] = arrayfun(@(tl,td)...
        block_turn_data(flyTracks.tStamps,flyTracks.rightTurns,tl,td),...
        t_light_start,[t_dark_start(2:end); flyTracks.tStamps(end)], 'UniformOutput', false);
end

% calculate turns metrics for both stimulus conditions
[light_seq,light_iti,light_switch,light_n] = ...
    cellfun(@(td,t) calculate_metrics(td,t,flyTracks.mazeOri),...
        light.turn_dat, light.t, 'UniformOutput', false);
[dark_seq,dark_iti,dark_switch,dark_n] = ...
    cellfun(@(td,t) calculate_metrics(td,t,flyTracks.mazeOri),...
        dark.turn_dat, dark.t, 'UniformOutput', false);

% convert data to arrays
light_n = cat(1,light_n{:});
light_switch = cat(2,light_switch{:})';
light_n(isnan(light_switch)) = 0;
dark_n = cat(1,dark_n{:});
dark_switch = cat(2,dark_switch{:})';
dark_n(isnan(dark_switch)) = 0;

% temp block masks
light_temp_block = cat(1,light.temp_block{:});
dark_temp_block = cat(1,dark.temp_block{:});

% calculate metrics for all condition blocks
clumpiness = NaN(8,flyTracks.nFlies);
switchiness = NaN(8,flyTracks.nFlies);
clumpiness_n = NaN(8,flyTracks.nFlies);
switchiness_n = NaN(8,flyTracks.nFlies);
right_bias = NaN(8,flyTracks.nFlies);
right_bias_n = NaN(8,flyTracks.nFlies);
for j=1:4
    % calculate switchiness for light blocks
    tmp_light_n = light_n(light_temp_block==j,:);
    block_weights = tmp_light_n./sum(tmp_light_n);
    switchiness(j*2-1,:) = nansum(light_switch(light_temp_block==j,:).*block_weights);
    switchiness(j*2-1,sum(tmp_light_n)<1) = NaN;
    switchiness_n(j*2-1,:) = sum(tmp_light_n);
    
    % calculate clumpiness for light blocks
    tmp_light_iti = light_iti(light_temp_block==j);
    tmp_light_iti = cat(2,tmp_light_iti{:});
    tmp_light_iti = cellfun(@(iti) cat(1,iti{:}), num2cell(tmp_light_iti,2), 'UniformOutput', false);
    tmp_light_n = cellfun(@numel,tmp_light_iti);
    clumpiness(j*2-1,:) = cellfun(@(iti) std(iti),tmp_light_iti)./(3600./tmp_light_n);
    clumpiness(j*2-1,sum(tmp_light_n)<1) = NaN;
    clumpiness_n(j*2-1,:) = tmp_light_n;
    
    % calculate right bias for light blocks
    tmp_light_seq = light_seq(light_temp_block==j);
    tmp_light_seq = cat(1,tmp_light_seq{:});
    tmp_light_seq = cellfun(@(iti) cat(1,iti{:}), num2cell(tmp_light_seq,1), 'UniformOutput', false);
    tmp_light_n = cellfun(@numel,tmp_light_seq);
    right_bias(j*2-1,:) = cellfun(@sum,tmp_light_seq)./tmp_light_n;
    right_bias(j*2-1,tmp_light_n<1) = NaN;
    right_bias_n(j*2-1,:) = tmp_light_n;
    
    % calculate switchiness for dark blocks
    tmp_dark_n = dark_n(dark_temp_block==j,:);
    block_weights = tmp_dark_n./sum(tmp_dark_n);
    switchiness(j*2,:) = nansum(dark_switch(dark_temp_block==j,:).*block_weights);
    switchiness(j*2,sum(tmp_dark_n)<1) = NaN;
    switchiness_n(j*2,:) = sum(tmp_dark_n);
    
    % calculate clumpiness for light blocks
    tmp_dark_iti = dark_iti(dark_temp_block==j);
    tmp_dark_iti = cat(2,tmp_dark_iti{:});
    tmp_dark_iti = cellfun(@(iti) cat(1,iti{:}), num2cell(tmp_dark_iti,2), 'UniformOutput', false);
    tmp_dark_n = cellfun(@numel,tmp_dark_iti);
    clumpiness(j*2,:) = cellfun(@(iti) std(iti),tmp_dark_iti)./(3600./tmp_dark_n);
    clumpiness(j*2,sum(tmp_dark_n)<1) = NaN;
    clumpiness_n(j*2,:) = tmp_dark_n;
    
    % calculate right bias for light blocks
    tmp_dark_seq = dark_seq(dark_temp_block==j);
    tmp_dark_seq = cat(1,tmp_dark_seq{:});
    tmp_dark_seq = cellfun(@(iti) cat(1,iti{:}), num2cell(tmp_dark_seq,1), 'UniformOutput', false);
    tmp_dark_n = cellfun(@numel,tmp_dark_seq);
    right_bias(j*2,:) = cellfun(@sum,tmp_dark_seq)./tmp_dark_n;
    right_bias(j*2,tmp_dark_n<1) = NaN;
    right_bias_n(j*2,:) = tmp_dark_n;
end

% assign outputs
out.right_bias.data = right_bias;
out.right_bias.n = right_bias_n;
out.switchiness.data = switchiness;
out.switchiness.n = switchiness_n;
out.clumpiness.data = clumpiness;
out.clumpiness.n = clumpiness_n;
out.labels = flyTracks.labels;
out.bouts.n = nBouts;
out.bouts.length = bout_length;
out.bouts.idx = bout_idx;
out.bouts.clumpiness = bout_clumpiness;
out.speed.data = speed_blocks;
out.speed.avg = mean_speed;


function [block_turns,block_t,block_temp] = block_turn_data(t_stamps,arm_dat,t1,t2)

mask = t_stamps >= t1 & t_stamps <= t2;
block_turns = arm_dat(mask,:);
block_t = t_stamps(mask);
if any(mask)
    block_temp = ceil(block_t(1)/3600);
else
    block_temp = NaN;
end


function [turn_sequence,iti,switchiness,n] = calculate_metrics(arm_dat,t,ori)

turn_idx = num2cell(~isnan(arm_dat),1);
arm_sequence = cellfun(@(aseq,ti) aseq(ti),num2cell(arm_dat,1), turn_idx, 'UniformOutput', false);
turn_sequence = cellfun(@score_turns,arm_sequence,num2cell(ori)','UniformOutput', false);
n = cellfun(@numel,turn_sequence);

% initialze placeholder data
switchiness = NaN(numel(n),1);
iti = cell(numel(n),1);

% Calculate turn metrics
no_turns = n < 2;
turn_t = cellfun(@(ti) t(ti), turn_idx(~no_turns), 'UniformOutput', false);
right_bias = cellfun(@sum,turn_sequence(~no_turns))./n(~no_turns);
switchiness(~no_turns) = cellfun(@(s,r,nt) sum((s(1:end-1)+s(2:end))==1)/(2*r*(1-r)*nt),...
    turn_sequence(~no_turns),num2cell(right_bias),num2cell(n(~no_turns)));
iti(~no_turns) = cellfun(@diff,turn_t,'UniformOutput',false);


function seq = score_turns(arm_seq,ori)

tSeq=diff(arm_seq);  
if ori
    seq = tSeq==1 | tSeq==-2;
else
    seq = tSeq==-1 | tSeq==2;
end


function speed = calculate_speed(cen,t)

dist = [NaN(1,size(cen,3));squeeze(sqrt(sum(diff(cen,1,1).^2,2)))];
speed = dist./[NaN;diff(t)];

