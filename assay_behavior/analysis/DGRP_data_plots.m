% Generates plots from the DGRP collection behavioral data.
% Load the DGRP data matrices (D_dgrp_behavior and D_dgrp_phys) before running.

% plot behavioral correlation matrix
plotCorr(D_dgrp_behavior.data,'Labels',D_dgrp_behavior.fields,'Patch',false);
axis('equal','tight');
title('DGRP behavioral corr mat');

% scree plot
subplot(1,2,2);
[~,~,nkeep] = plot_pca_bootstrap(D_dgrp_behavior.data,100,95,'noncummulative',[.7 .7 .7]);
set(gca,'YLim',[1E-2 1E2],'XLim',[1 size(D_dgrp_behavior.data,2)]);
title(sprintf('thresh = %i',nkeep));

%% plots for physiological/molecular measurements

% plot behavioral correlation matrix
plotCorr(D_dgrp_phys.data,'Labels',D_dgrp_phys.fields,'Patch',false);
axis('equal','tight');
title('DGRP physiological corr mat');

% scree plot
subplot(1,2,2);
[~,~,nkeep] = plot_pca_bootstrap(D_dgrp_phys.data,100,95,'noncummulative',[.7 .7 .7]);
set(gca,'YLim',[1E-2 1E2],'XLim',[1 numel(D_dgrp_phys.fields)]);
title(sprintf('thresh = %i',nkeep));

%% combined analysis on both behavioral and physiological data

% dimensionality plots (connected components log-hist)
figure;
D = [D_dgrp_behavior; D_dgrp_phys];
d_label = {'DGRP behavior data';'DGRP physiological data'};
lhs = gobjects(numel(D),1);
for i=1:numel(D)
    subplot(1,2,i);
    d = D(i).data;
    nfeat = size(d,2);
    nreps = 100;
    nflies = size(d,1);
    
    comp_cts = zeros(nreps);
    for j=1:nreps
        if mod(j,10)==0
           fprintf('Conn. comp. simulation %i of %i\n',j,nreps); 
        end
        dd = d(randi(nflies,[nflies 1]),:);

        % compute conn comp spectrum
        comp_cts(j,:) = conn_comp_spectrum(corr(dd),nreps);
    end
    
    cts = histc(comp_cts(:),1:nfeat);
    lhs(i) = bar(1:nfeat,log10(cts),'FaceColor',[.65 .65 .65],'EdgeColor','none');
    xlabel('conn. comp');
    ylabel('log(counts)');
    drawnow;
    ah = gca;
    ah.Units = 'inches';
    ah.Position(3:4) = [2 1];
    title(d_label{i});
    drawnow;
end

%%

% DGRP t-SNE plots
fh = figure;
dist_type = 'euclidean';
perplex = 8;
batch = {'D12';'D3'};
leg_loc = 'NorthEast';
fdir = 'D:\decathlon_data_and_analysis\decathlon_analysis\figures\multidimensional_scaling\';
fname = sprintf('BABAM_tsne_perplex%i_%s.fig',perplex,dist_type);
ncol = 2;
opts = statset;
opts.MaxIter = 5000;
opts.TolFun = 1E-15;
buffer_scale = 1;
alg = 'exact';

ah = subplot(1,ncol,1);

% caluclate fraction missing
z_data = zscore(D_dgrp_phys.data);
embedded_data = tsne(z_data,'Distance',dist_type,...
    'Perplexity',perplex,'Algorithm',alg,'Options',opts);
pretty_scatter(embedded_data(:,1),embedded_data(:,2),'k');
title(sprintf('t-SNE %s - lines physiological',dist_type));
axis('equal');
lim = std(embedded_data(:))*3.5;
set(gca,'XLim',[-lim lim].*buffer_scale,'YLim',[-lim lim].*buffer_scale);

% concatenate grouped data
ah = subplot(1,ncol,2);
embedded_data = tsne(z_data','Distance',dist_type,...
    'Perplexity',perplex,'Algorithm',alg,'Options',opts);
pretty_scatter(embedded_data(:,1),embedded_data(:,2),'k');
title(sprintf('t-SNE %s - metrics physiological',dist_type));
axis('equal');
lim = std(embedded_data(:))*3.5;

set(gca,'XLim',[-lim lim].*buffer_scale,'YLim',[-lim lim].*buffer_scale);

% caluclate fraction missing
figure;
z_data = zscore(D_dgrp_behavior.data);
ah = subplot(1,ncol,1);
embedded_data = tsne(z_data,'Distance',dist_type,...
    'Perplexity',perplex,'Algorithm',alg,'Options',opts);
pretty_scatter(embedded_data(:,1),embedded_data(:,2),'k');
title(sprintf('t-SNE %s - lines behavior',dist_type));
axis('equal');
lim = std(embedded_data(:))*3.5;

set(gca,'XLim',[-lim lim].*buffer_scale,'YLim',[-lim lim].*buffer_scale);

% concatenate grouped data
ah = subplot(1,ncol,2);
embedded_data = tsne(z_data,'Distance',dist_type,...
    'Perplexity',perplex,'Algorithm',alg,'Options',opts);
pretty_scatter(embedded_data(:,1),embedded_data(:,2),'k');
title(sprintf('t-SNE %s - metrics behavior',dist_type));
axis('equal');
lim = std(embedded_data(:))*3.5;

set(gca,'XLim',[-lim lim].*buffer_scale,'YLim',[-lim lim].*buffer_scale);
