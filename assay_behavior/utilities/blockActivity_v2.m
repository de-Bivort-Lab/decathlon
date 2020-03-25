function [varargout] = blockActivity_v2(spd)

% blockActivity divides margo speed traces into discreet bouts
%
% input:
%   s_map           ->  memmap to raw speed data file (ie. expmt.data.speed.map)
%
% outputs:
%
%   block_indices   ->  nBouts x 2 cell array of frame indices where bout
%                       transitions occurred
%   lag_thresh      ->  autocorrelation lag threshold value in number of
%                       frames for defining duration of movement bouts
%   speed_thresh    ->  threshold value acquired by fitting two-component 
%                       gmm to log(speed) and computing intersection
%
% compute autocorrelation and find conservative 
% cutoff for bout discretization
% 

if size(spd,1) < 50000
    smpl = 1:size(spd,1);
else
    smpl = 1:50000;
end

% find speed threshold cutoff from log speed
[intersect,class_means] = dec_kthresh_distribution(log(medfilt1(spd(:),3)));
speed_thresh = exp(intersect);
s = spd;
s(isnan(s)) = 0;

% plot results
fh = figure('Name','Bout Parsing');
ah1 = subplot(1,2,1);
hold on;
x = linspace(-6,6,60);
ss = log(medfilt1(spd(:),3,[],1));
ss(ss<log(1E-4)) = [];
y = ksdensity(ss(randperm(numel(ss),10000)),x);
plot(x,y);
max_y = max(y);
max_y(isnan(max_y)|max_y<0.01)=0.01;
plot([intersect intersect],[0 max_y*1.2],'r--');
plot([class_means(1) class_means(1)],[0 max_y*1.2],'k--');
plot([class_means(2) class_means(2)],[0 max_y*1.2],'k--');
set(ah1,'XLim',[-6 6],'YLim',[0 max_y*1.25]);
xlabel('log(speed)');
title('speed threshold');

% find frames where transitioned from 
moving = s > speed_thresh;
f_moving = medfilt1(single(moving),3,[],1);
f_moving = f_moving(smpl,:);

% compute autocorrelation
ac = acf(f_moving(:),100);
y = medfilt1(ac,3);
dy = [NaN; diff(y)];
[~,lag_thresh] = min(dy);

% plot results
ah2 = subplot(1,2,2);
plot(y,'k-','LineWidth',1.5);
hold on; plot([lag_thresh lag_thresh],[0 1],'r--');
xlabel('lag (frames)');
title('speed autocorrelation threshold');

% median filter data by lag_thresh/2 to discretize bouts
if (lag_thresh>6) 
    s = medfilt1(s,round(lag_thresh/2),[],1);
else
    s = medfilt1(s,3,[],1);
end

moving = s > speed_thresh;
f_moving = int8(medfilt1(single(moving),3,[],1));
transitions = diff(f_moving,1,1);
clear moving
transitions = cat(1,int8(zeros(1,size(transitions,2))),transitions);
transitions = num2cell(transitions,1);

% get activity bout stops and starts
stops = cellfun(@(x,y) find(x==-1)-1, transitions,'UniformOutput',false);
starts = cellfun(@(x,y) find(x==1), transitions,'UniformOutput',false);

% free speed map
clear s moving transitions

% filter by bout lengths
block_indices = arrayfun(@(x,y) ...
                filterShortBouts(x,y,lag_thresh/2),...
                starts,stops,'UniformOutput',false);
            
            
for i=1:nargout
    switch i
        case 1, varargout{i} = block_indices;
        case 2, varargout{i} = lag_thresh;
        case 3, varargout{i} = speed_thresh;
        case 4, varargout{i} = fh;
    end
end



    
function idx = filterShortBouts(starts,stops,duration)

    % discard last bout if it starts but doesn't end
    if any(size(starts{:}) ~= size(stops{:}))

        start_end = length(starts{:}) > length(stops{:});
       
        if start_end
            starts{:}(end) = [];
        else
            stops{:}(end) = [];
        end
    end
    

    if ~isempty(starts) && ~isempty(stops)
        bout_length = abs(starts{:}-stops{:});
        long_bout = bout_length >= duration;
        idx = [starts{:}(long_bout) stops{:}(long_bout)];

        if size(idx,1)>1
            % ensure that lower index comes first
            [~,i]=min(idx,[],2);
            idx(i==2,:)=[];
        end
    else
        idx = [];
    end
    
    idx = uint32(idx);
    
    
    

        
        