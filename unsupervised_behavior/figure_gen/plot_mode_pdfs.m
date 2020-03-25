function plot_mode_pdfs(idxMap,pdfs)

figure;
ah1 = gca;
[~,p] = sort(sum(pdfs));
p = fliplr(p);
[~,pp] = sort(median(pdfs,2));
if size(pdfs,1)>1
    z = linkage(pdfs,'average','correlation');
    [~,~,zp] = dendrogram(z,0);
else
    zp = 1;
end

sorted_pdfs = pdfs(zp,p);
imagesc(sorted_pdfs');
colors = [1 1 1; .75 .75 .75; .5 0 .8; 1 .1 0; 1 .9 0];
cmap=interp1([1 51 102 153 205],colors,1:205);
colors = [1 1 1; 0 0 0];
cmap=interp1([0 255],colors,0:255);
cmap = interp1([0 102 205],[0 0 0; 1 0 0; 1 1 0],0:205);
cmap = bone;
cmap = bone;
cmap(size(cmap,1):-1:1,:) = cmap;
colors = [255 255 255; 149 174 181; 64 64 89; 0 0 0]./255;
cmap = interp1([0 45 128 205],colors,0:205);



colorbar
colormap(ah1,cmap);
caxis(ah1,[0 prctile(pdfs(:),99)]);
title('behavioral mode PDFs');
xlabel('individual flies');
ylabel('mode no.');

% get unique modes
unique_modes = unique(idxMap);
unique_modes(unique_modes==0)=[];

% count a fly as visiting a mode at least once if move than 10 frames are
% labeled as a particular identity
figure;
total_pdf = sum(pdfs,1)./sum(pdfs(:));
fpm = zeros(size(idxMap));
modeCentroid = NaN(numel(unique_modes),2);
idxMap_bordered = idxMap;
pts = 256;
% cmap = interp1([0 pts*(1/3) pts*(2/3) pts],...
%     [1 1 1; .1 .1 1; 1 1 0; 1 0 0],0:pts);
% cmap = interp1([0 pts*(1/3) pts*(2/3) pts-1E-8 pts],...
%     [1 1 1; .1 .1 1; 1 1 0; 1 0 0; 0 0 0],0:pts);
colors = [1 1 1; 0 0 0; 0 0 0];
cmap=interp1([1 204 205],colors,1:205);
cmap = interp1([0 102 204 205],[0 0 0; 1 0 0; 1 1 0; 1 1 1],0:205);
cmap = bone;
colors = [255 255 255; 149 174 181; 64 64 89; 0 0 0]./255;
cmap = interp1([0 45 128 205],colors,0:205);


for i=1:numel(unique_modes)
    modeMask = idxMap_bordered==unique_modes(i);
    props = regionprops(modeMask,'Centroid');
    if ~isempty(props)
        modeCentroid(i,:) = props(1).Centroid;
    end
    fpm(modeMask)=total_pdf(i);
end

molaspass=interp1([1 51 102 153 204 205],[0 0 0; 0 0 .75; .5 0 .8; 1 .1 0; 1 .9 0; 1 1 1],1:205);
y_max = max(total_pdf)*(205/204);

% ensure boundaries are set to blank
fpm(idxMap==0)=y_max;
idx_map_border = imerode(idxMap<1,strel('disk',2));
fpm(idx_map_border) = 0;
imagesc(fpm);
axis equal tight off
colormap(cmap);
colorbar
% for i=1:length(modeCentroid)
%     text(modeCentroid(i,1),modeCentroid(i,2),...
%         num2str(unique_modes(i)),'Color','k',...
%         'HorizontalAlignment','center','FontSize',8);
% end
xlabel('PDF');
set(gca,'XLim',[18 497],'YLim',[47 497]);