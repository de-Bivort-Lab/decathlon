function plot_watershed_map(idx_map,density,varargin)

numbered = false;
ah = gca;
for i=1:numel(varargin)
    arg = varargin{i};
    if ischar(arg)
        switch arg
            case 'Parent'
                i=i+1;
                ah = varargin{i};
            case 'CLim'
                i=i+1;
                clim = varargin{i};
            case 'Numbered'
                i=i+1;
                numbered = varargin{i};
        end
    end
end

mask = density < 1E-6;
big_mask = imdilate(~mask,strel('disk',2));
plot_map = idx_map > 0;
plot_map(~big_mask) = true;
imagesc(plot_map,'Parent',ah);
hold(ah,'on');
axis(ah,'equal','tight','off');
colormap(ah,'gray');

if numbered
   unique_modes = unique(idx_map(:));
   unique_modes = 1:unique_modes(end);
   for i=1:numel(unique_modes)
      idx_mask = idx_map == unique_modes(i);
      cen = regionprops(idx_mask,density,'WeightedCentroid');
      cen = cen.WeightedCentroid;
      text(cen(1),cen(2),sprintf('%i',i),'Color','k','FontSize',10,...
          'HorizontalAlignment','center','HitTest','off','PickableParts','none');
   end
end