function speed = average_3frame_speed(cen,ifi)

% calculate a 3 frame average of position
cen_previous = [NaN(1,2,size(cen,3)); cen(1:size(cen,1)-1,:,:)];
cen_current = cen;
cen_next = [cen(2:size(cen,1),:,:); NaN(1,2,size(cen,3))];
cen = cat(4,cen_previous,cen_current,cen_next);
cen = nanmean(cen,4);

% calculate speed from this average
dist = [NaN(1,size(cen,3));squeeze(sqrt(sum(diff(cen,1,1).^2,2)))];
speed = dist./ifi;
