clearvars -except RNAseq

% load decathlon data
D = load_decathlon_structs('D13_als');

% generate distilled matrix and restrict to BerlinK-iso decathlon-2
opts = {'CollapseMode','PCA','CollapseFields','all',...
    'Standardize',true,'ImputeMode','none'};
D = pair_decathlon_structs(D,opts{:});
D = standardize_by_field(D);
D = D(1);
d = D.data(191:end,:);

% quantile normalize the reads
norm_reads = quantile_normalize(RNAseq.reads);


%% plot correlation matrices

% individual flies
ahs = subplot_array(6);
plot_num = 1;
title(ahs(plot_num),'flies x flies behavior correlation matrix');
plotCorr(d','Parent',ahs([plot_num plot_num+3]),'Patch',false);
axis(ahs(plot_num),'off');
delete(ahs(plot_num+3));
ahs(plot_num+3) = subplot(2,3,6);
ahs(plot_num+3).YLim(2) = 100;
[~,~,k] = plot_pca_bootstrap(d',100,95,'noncummulative',[1 0 0]);
ahs(plot_num+3).YLim = [0.01 100];
title(sprintf('individual flies, k = %i',k));

% filter out data from empty individuals
is_sequenced = ~all(isnan(norm_reads));
norm_reads = norm_reads(:,is_sequenced);
d = d(is_sequenced,:);

% flies by flies gene expression
plot_num = 2;
title(ahs(plot_num),'flies x flies behavior correlation matrix');
plotCorr(norm_reads,'Parent',ahs([plot_num plot_num+3]),'Patch',false);
axis(ahs(plot_num),'off');
delete(ahs(plot_num+3));
ahs(plot_num+3) = subplot(2,3,plot_num+3);
[~,~,k] = plot_pca_bootstrap(norm_reads,100,95,'noncummulative',[1 0 0]);
ahs(plot_num+3).YLim = [0.01 100];
title(sprintf('individual flies, k = %i',k));

% genes by genes
plot_num = 3;
title(ahs(plot_num),'flies x flies behavior correlation matrix');
r = corr(norm_reads,'Type','Spearman');
r(isnan(r))=0;

% sort rows and columns by hierarchical clustering
Z=linkage(r);
f=figure;
[ZH, ZT, Zoutperm]=dendrogram(Z,0);
close(f);
r=r(Zoutperm,Zoutperm);


axis(ahs(plot_num),'off');
delete(ahs(plot_num+3));
ahs(plot_num+3) = subplot(2,3,plot_num+3);
ahs(plot_num+3).YLim(2) = 100;
[~,~,k] = plot_pca_bootstrap(norm_reads',100,95,'noncummulative',[1 0 0]);
ahs(plot_num+3).YLim = [0.01 100];
title(sprintf('individual flies, k = %i',k));

%% Lasso regression model

[~,~,k] = plot_pca_bootstrap(norm_reads(~all(isnan(norm_reads),2),:)',20,95,'noncummulative',[1 0 0]);
set(gca,'YLim',[0.01 100]);
[~,Y] = pca(d);
X = norm_reads';
X(:,all(X==0));

out = cell(size(Y,2),1);
for i=1:numel(out)
    out{i}=decathlonLasso(X,Y(:,i));
end

%% Gene expression NxD matrix plot

r = corr(norm_reads,'Type','Spearman');
r(isnan(r))=0;
Z=linkage(r);
[ZH, ZT, Zoutperm]=dendrogram(Z,0);

im = norm_reads(:,Zoutperm);
[~,perm] = sort(nanmean(im,2));
figure; imagesc(log(im(perm,:)'));
colorbar
colormap(molaspass)
axis('tight')
ylabel('individual flies');
xlabel('genes');
set(gca,'XTick',[],'YTick',[])