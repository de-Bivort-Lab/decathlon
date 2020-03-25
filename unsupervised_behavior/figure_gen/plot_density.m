function plot_density(density,idx_map,varargin)

ah = [];
clim = [0 max(density(:))];
numbered = false;
outline_density = density;
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
            case 'OutlineDensity'
                i=i+1;
                outline_density = varargin{i};
        end
    end
end

if isempty(ah)
    ah = gca;
end

% plot the density
imagesc(density,'Parent',ah);
axis(ah,'equal','tight','off');
pts = 255;
cmap = interp1([0 pts*(1/3) pts*(2/3) pts-1E-8 pts],...
    [1 1 1; .1 .1 1; 1 1 0; 1 0 0; 0 0 0],0:pts);
colors = [1 1 1; 0 0 .75; .5 0 .8; 1 .1 0; 1 .9 0; 0 0 0];
cmap=interp1([1 51 102 153 204 205],colors,1:205);
cmap = interp1([0 51 204 205],[1 1 1; .5 .5 .5; 0 0 0; .6 0 .6],0:205);
cmap = interp1([0 102 204 205],[0 0 0; 1 0 0; 1 1 0; 1 1 1],0:205);
colors = [0 0 0; 140 77 134; 203 162 179; 255 255 255]./255;
cmap = interp1([0 115 170 205],colors,0:205);
cmap(size(cmap,1):-1:1,:) = cmap;
cmap = bone;
cmap(size(cmap,1):-1:1,:) = cmap;
colors = [255 255 255; 149 174 181; 64 64 89; 0 0 0]./255;
cmap = interp1([0 45 128 205],colors,0:205);

caxis(clim);
colormap(ah,cmap);
colorbar(ah);

% plot the watershed mask
mask = outline_density<1E-6;
mask_edge = edge(mask,'canny');
binim = idx_map == 0;
binim = binim & ~mask;
mode_outlines = binim | mask_edge;
idx_map(mask) = 0;

hold(ah,'on');
imh2 = imagesc(imdilate(mode_outlines,strel('disk',1)),'Parent',ah);
imh2.AlphaData = imdilate(mode_outlines,strel('disk',1));

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
set(ah,'XLim',[18 497],'YLim',[47 497]);