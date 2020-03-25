function [starts,stops,durations,sampleFrames,sampleDurations] = modeBouts(fID,modes,tDur)

% create mask for each mode in fID
modes = (1:modes(end))';
[starts,stops,durations,sampleFrames,sampleDurations]  = ...
    arrayfun(@(x) getBouts(x,fID,tDur),modes,'UniformOutput',false);



function [starts,stops,durations,sampleFrames,sampleDurations] = getBouts(targetMode,fID,tDur)

istarget = fID==targetMode;
transitions = [0;diff(istarget)];
starts = find(transitions == 1);
stops = find(transitions == -1);

% discard incomplete bouts
if ~isempty(starts) && ~isempty(stops)...
        && stops(1) < starts(1)
    stops(1) = [];
end
if ~isempty(starts) && ~isempty(stops)...
        && starts(end) > stops(end)
    starts(end) = [];
end

% get bout durations
durations = stops - starts;
sampleFrames=starts(durations>tDur);
sampleDurations = durations(durations>tDur);