function stable_mode_mask = filter_speed_mask(speed_mask,min_length)

trans = [0;int8(diff(speed_mask))];
starts = find(trans==1);
stops = find(trans==-1);

if stops(1)<starts(1)
    stops(1) = [];
end
if starts(end) > stops(end)
    starts(end) = [];
end

bout_length = stops - starts;
short_bouts = bout_length < min_length;
starts(short_bouts) = [];
stops(short_bouts) = [];

bout_idx = arrayfun(@(sa,sp) sa:sp, starts, stops, 'UniformOutput', false);
bout_idx = cat(2,bout_idx{:})';

stable_mode_mask = false(size(speed_mask));
stable_mode_mask(bout_idx) = true;