% Generates plots from the BABAM behavioral data (Robie et. al, Cell 2015)
% Load the BABAM data matrix (D_babam) before running.

data = D_babam.data;

figure;
plotCorr(data,'Patch',false);
axis('equal','tight');
title('BABAM corr mat');

% plot cross validation results
k_folds = 1;
figure;
subplot(1,2,1);
[tr,te,~] = cross_validate_pca(data,'KFolds',k_folds,'TestSize',0.1);
[~,mi] = min(mean(te,2));
title(sprintf('thresh = %i',mi)); 

% scree plot
subplot(1,2,2);
[~,~,nkeep] = plot_pca_bootstrap(z_data,20,95,'noncummulative',[.7 .7 .7]);
set(gca,'YLim',[1E-2 1E2]);
title(sprintf('thresh = %i',nkeep));

% compute the number of effective dimensions via N-drop scree intersection
max_reps = 150;
z_data = zscore(data);

plot_pca_n_drop_scree(z_data,max_reps);
title({'Robie Data';'PCA drop N 150 reps'});

% compute effective dimensionality connected components in correlation matrix
figure; hold on;
d = z_data;
[out,intrinsic_dim,~]=decathlonConnCompSweep(d,1000);
bins = 1:8:size(d,2);
cts = histc(intrinsic_dim,bins);
bar(bins,log(cts),'EdgeColor','none','FaceColor',[0 0 0]);
ylabel('log(count)');
xlabel('effective dimensionality');
set(gca,'XLim',[0 size(d,2)]);


% BABAM t-SNE plots
fh = figure;
dist_type = 'euclidean';
perplex = 20;
fdir = 'D:\decathlon_data_and_analysis\decathlon_analysis\figures\multidimensional_scaling\';
fname = sprintf('BABAM_tsne_perplex%i_%s.fig',perplex,dist_type);
ncol = 2;
opts = statset;
opts.MaxIter = 5000;
opts.TolFun = 1E-12;
opts.Display = 'iter';

ah = subplot(1,ncol,1);

% tSNE plot for individuals
embedded_data = tsne(data,'Distance',dist_type,'Perplexity',perplex,'Options',opts);
pretty_scatter(embedded_data(:,1),embedded_data(:,2),'k');
title({'Robie data';sprintf('t-SNE %s - lines',dist_type)});
axis('equal');
set(gca,'XLim',[-70 70],'YLim',[-70 70]);

% tSNE plot for behavioral metrics
ah = subplot(1,ncol,2);
embedded_data = tsne(data','Distance',dist_type,'Perplexity',perplex,'Algorithm','exact','Options',opts);
pretty_scatter(embedded_data(:,1),embedded_data(:,2),'k');
title({'Robie data';sprintf('t-SNE %s - metrics',dist_type)});
axis('equal');
set(gca,'XLim',[-80 80],'YLim',[-80 80]);

