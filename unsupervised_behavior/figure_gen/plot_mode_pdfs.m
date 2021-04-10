function plot_mode_pdfs(idxMap,pdfs,mask)

figure;
for i=1:2
    ah = subplot(2,1,i);
    if mod(i,2)
        d = pdfs(mask,:);
        title_str = 'inbred';
    else
        d = pdfs(~mask,:);
        title_str = 'outbred';
    end
    [~,p] = sort(sum(d));
    p = fliplr(p);
    [~,pp] = sort(median(d,2));
    if size(d,1)>1
        z = linkage(d,'average','correlation');
        [~,~,zp] = dendrogram(z,0);
    else
        zp = 1;
    end

    sorted_pdfs = d(zp,p);
    imagesc(sorted_pdfs');
    colors = [255 255 255; 149 174 181; 64 64 89; 0 0 0]./255;
    cmap = interp1([0 45 128 205],colors,0:205);

    colorbar
    colormap(ah,cmap);
    caxis(ah,[0 prctile(d(:),99)]);
    title(sprintf('behavioral mode PDFs (%s)',title_str));
    xlabel('individual flies');
    ylabel('mode no.');
end

% get unique modes
unique_modes = unique(idxMap);
unique_modes(unique_modes==0)=[];

% count a fly as visiting a mode at least once if move than 10 frames are
% labeled as a particular identity
figure;

for j=1:2
    ah = subplot(2,1,j);
    if mod(j,2)
        d = pdfs(mask,:);
        title_str = 'inbred';
    else
        d = pdfs(~mask,:);
        title_str = 'outbred';
    end
    total_pdf = sum(d,1)./sum(d(:));
    fpm = zeros(size(idxMap));
    modeCentroid = NaN(numel(unique_modes),2);
    idxMap_bordered = idxMap;
    pts = 256;
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

    y_max = max(total_pdf)*(205/204);

    % ensure boundaries are set to blank
    fpm(idxMap==0)=y_max;
    idx_map_border = imerode(idxMap<1,strel('disk',2));
    fpm(idx_map_border) = 0;
    imagesc(fpm);
    axis equal tight off
    colormap(cmap);
    colorbar
    xlabel('PDF');
    set(gca,'XLim',[18 497],'YLim',[47 497]);
    title(sprintf('mode map (%s)',title_str));
end