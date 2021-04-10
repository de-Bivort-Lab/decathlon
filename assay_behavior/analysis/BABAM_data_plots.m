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

%% compute effective dimensionality connected components in correlation matrix

figure;
d = z_data;
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
bar(1:nfeat,log10(cts),'FaceColor',[.65 .65 .65],'EdgeColor','none');
xlabel('conn. comp');
ylabel('log(counts)');
drawnow;
ah = gca;
ah.Units = 'inches';
ah.Position(3:4) = [2 1];
drawnow;



%% BABAM t-SNE plots
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

