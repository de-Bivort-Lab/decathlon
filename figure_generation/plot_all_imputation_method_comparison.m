
fpath = ['D:\decathlon_data_and_analysis\decathlon_analysis\matrices\' ...
    'decathlon_paper\nuisance unregressed\D123_olfaction_timeofday_added.mat'];

als_data = cell(4,1);
ncols = ceil(sqrt(numel(als_data)));
nrows = ceil(numel(als_data)/ncols);
bfh = figure;
for i=1:numel(als_data)
    load(fpath);
    D = impute_decathlon_structs(D,'ImputeMode','als');
    D = standardize_by_field(D);
    als_data{i} = D;
    [fh1, fh2] = plot_paired_corrmats(D,'CollapseFields','circadian');
    close(fh1);
%     figure(bfh);
%     subplot(nrows,ncols,i);
%     bootstrap_decathlon_structs(D,50);
    drawnow
    pause(0.5);
    clear D
end

%%

ncols = numel(als_data);
cmap=interp1([1 47 128 129 169 256],...
    [0 1 1; 0 .2 1; 0 0 0; 0 0 0 ; 1 .1 0; 1 1 0],1:256);
d = pair_decathlon_structs(als_data{i},'CollapseFields','Circadian');
[p,~] = get_cluster_permutation(d);
labels = cell(numel(p),1);
labels(:) = {''};
pf = figure;
pa = gca;
cf = figure;

for i=1:numel(als_data)
    d = pair_decathlon_structs(als_data{i},'CollapseFields','Circadian');
    d = apply_cluster_permutation(p, d);
    for j=1:3
        figure(cf);
        ah = subplot(3,ncols,(j-1)*ncols + i);
        plotCorr(d(j).data,'Labels',labels,'Cluster',false,'Parent',[ah pa],'Patch',false);
        title(sprintf('D%i - als filled',j));
        axis equal tight off
    end
end
close(pf);


%%

opts = {'Marker'; 'o'; 'LineStyle'; 'none'; 'MarkerEdgeColor'; 'none';...
        'MarkerSize'; 1; 'LineWidth'; 2};
pairs = unique_idx_pairs(numel(als_data),1);
ncols = size(pairs,1);
%colors = rand(size(pairs,1),3);

figure;

for j=1:3
    for i=1:ncols
        ah = subplot(3,size(pairs,1),(j-1)*ncols + i);
        d1 = als_data{pairs(i,1)};
        d2 = als_data{pairs(i,2)};
        f = d1(j).imputed;
        %plot(d1(j).data(f),d2(j).data(f),'MarkerFaceColor','k',opts{:});
        pts = linspace(-6,6,100);
        [x,y] = meshgrid(pts,pts);
        kde = ksdensity([d1(j).data(f),d2(j).data(f)],[x(:) y(:)]);
        imagesc(reshape(kde,numel(pts),numel(pts)))
        axis tight equal
        %set(ah,'XLim',[-6 6],'YLim',[-6 6]);
        ticks = linspace(1,numel(pts),7);
        set(ah,'XTick',ticks,'XTickLabels',-6:2:6,'YTick',ticks,'YTickLabels',-6:2:6);
        title(sprintf('D%i',j));
        xlabel(sprintf('run #%i',pairs(i,1)));
        ylabel(sprintf('run #%i',pairs(i,2)));
    end

end

%%

D_avg = D;
for i=1:numel(D)
    d_mats = cellfun(@(ad) ad(i).data, als_data, 'UniformOutput', false);
    d_mats = cat(3,d_mats{:});
    D_avg(i).data = mean(d_mats,3);
end
