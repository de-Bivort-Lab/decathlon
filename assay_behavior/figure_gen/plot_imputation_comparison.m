% Compare imputation methods with simulated ground truth data.
% Used to generate panels of supplemental fig. 4

% define data dimensions
m = 400;        % number of observations
n = 20;        % number of features
q = 3;         % ground truth dimensionality
p = 5;          % number of features in each correlated cluster
cluster_cov = 0.5;

% impart correlations to cluster by adding shared noise
sim_cov = zeros(n);
for i=1:q
    j = (i-1)*p;
    if i==q
        sim_cov(j+1:j+p,j+1:j+p) = cluster_cov/2;
    else
        sim_cov(j+1:j+p,j+1:j+p) = cluster_cov;
    end
end
sim_cov(logical(diag(ones(n,1))))=1;
sim_data = mvnrnd(zeros(n,1),sim_cov,m);
sim_data = zscore(sim_data);

% display covariance matrix
figure;
ah=subplot(1,2,1);
imagesc(cov(sim_data));
k = n-q*(p-1);
title({sprintf('num features = %i, clusters = %i',n,q);...
    sprintf('effective dimensions = %i',k)});
colormap(ah,egoalley);
colorbar;

% plot connected components output
out = decathlonConnCompDropN(sim_data);
ah=subplot(1,2,2);
imagesc(out);
hold on;
plot([k k],[0.5 size(out,1)+.5],'k--');
ylabel('features dropped');
xlabel('effective dimensionality');
title('connected components heatmap');
colormap(ah,logjet_cmap);

%% test different methods of imputation

g_truth = nanzscore(sim_data);
frac_missing = linspace(0.1,0.5,9);
nreps = 100;

% mean imputation
mean_impute_mse = zeros(nreps,numel(frac_missing));
for i=1:numel(frac_missing)
    
    fprintf('iteration %i of %i\n',i,numel(frac_missing));
    for j=1:nreps
       D.data = g_truth;
       mask = rand(size(g_truth))<frac_missing(i);
       D.data(mask) = NaN;
       D.data = nanzscore(D.data);
       opts = {'Standardize',false};
       D = impute_decathlon_structs(D,opts{:});

       % calculate error
       residuals = D.data(mask) - g_truth(mask);
       mean_impute_mse(j,i) = mean(residuals.^2);
    end
end

% ALS median imputation
g_truth = nanzscore(sim_data);
frac_missing = linspace(0.1,0.5,9);
nreps = 1;

% mean imputation
mse = zeros(nreps,numel(frac_missing));
data = cell(nreps,numel(frac_missing));
avg_data = cell(nreps,numel(frac_missing));
masks = data;
als_data = data;
avg_mse = mse;
for i=1:numel(frac_missing)
    
    fprintf('iteration %i of %i\n',i,numel(frac_missing));
    for j=1:nreps
       D.data = g_truth;
       mask = rand(size(g_truth))<frac_missing(i);
       masks{j,i} = mask;
       D.data(mask) = NaN;
       [D,als_data{i,j}] = avg_als_impute(D,200);
       data{j,i} = D.data;

       % calculate error
       residuals = D.data(mask) - g_truth(mask);
       mse(j,i) = mean(residuals.^2);
       
       als_avg = cellfun(@(d) mean(d,3), als_data{i}, 'UniformOutput',false);
       avg_mse(i) = mean((als_avg{1}(masks{i}) - g_truth(masks{i})).^2);
       avg_data(i) = als_avg;
    end
end

% plot MSE between ALS/mean imputed data and the ground truth
figure;
clust_feat = false(size(mask));
clust_feat(:,1:n <= p*q) = true;
clust_mse = NaN(numel(frac_missing),2);
for i=1:numel(frac_missing)
    clust_mse(i,1) = mean((data{i}(masks{i}&clust_feat) - g_truth(masks{i}&clust_feat)).^2);
    clust_mse(i,2) = mean((data{i}(masks{i}&~clust_feat) - g_truth(masks{i}&~clust_feat)).^2);
end
hold on;
lh = gobjects(3,1);
lh(1)=plot(frac_missing,clust_mse(:,1),'r-');
lh(2)=plot(frac_missing,clust_mse(:,2),'k');
lh(3)=plot(frac_missing,mean(mean_impute_mse,1),'k--');
ylabel('MSE');
xlabel('fraction missing');
title('imputation error');
legend(lh,{'cluster ALS';'non-cluster ALS';'mean imputed'},'Location','SouthEast');
set(gca,'YLim',[0 1.2]);

% scatter ALS imputed data vs ground truth
figure;
for i=1:numel(frac_missing)
    
    subplot(3,3,i);
    
    % cluster features
    infilled = data{i}(masks{i}&clust_feat);
    gt = g_truth(masks{i}&clust_feat);
    pretty_scatter(gt,infilled,'r','MarkerSize',2); hold on;
    r1 = corr([gt infilled],'Type','pearson');
    
    % non cluster features
    infilled = data{i}(masks{i}&~clust_feat);
    gt = g_truth(masks{i}&~clust_feat);
    r2 = corr([gt infilled],'Type','pearson');
    pretty_scatter(gt,infilled,'k','MarkerSize',2);
    ylabel('imputed');
    xlabel('ground truth');
    title(sprintf('frac. missing = %0.2f, r = %0.2f',frac_missing(i),r1(1,2)));
    axis('equal');
    legend({sprintf('cluster r=%0.2f',r1(1,2));...
        sprintf('non-cluster r=%0.2f',r2(1,2))});
    set(gca,'XLim',[-3.5 3.5],'YLim',[-3.5 3.5]);
end


% generate scree plots for ground truth, ALS, and mean imputed data
figure;
for i=1:numel(frac_missing)
    subplot(3,3,i);
    d = data{1,i};
    mask = rand(size(g_truth))<frac_missing(i);
    D.data = g_truth;
    D.data(mask) = NaN;
    opts = {'Standardize',false};
    D = impute_decathlon_structs(D,opts{:});
    [~,~,~] = plot_pca_bootstrap(d,100,95,'noncummulative',[1 0 0]);
    [~,~,~] = plot_pca_bootstrap(g_truth,100,95,'noncummulative',[0 0 1]);
    [~,~,~] = plot_pca_bootstrap(D.data,100,95,'noncummulative',[0 1 1]);
end

%% simulate several conn comp plots

% define data dimensions
m = 500;        % number of observations
n = 30;        % number of features
q = [3 6];         % ground truth dimensionality
p = [2 5];          % number of features in each correlated cluster
cluster_cov = [0.25 0.5];
fhs = gobjects(3,1);
for i=1:numel(fhs)
   fhs(i) = figure; 
end

for ii=1:numel(q)
    for jj=1:numel(p)
        for kk=1:numel(cluster_cov)
            
            ct = (ii-1)*numel(q)*numel(p) + (jj-1)*numel(p) + kk;
            
            % impart correlations to cluster by adding shared noise
            sim_cov = zeros(n);
            for i=1:q(ii)
                j = (i-1)*p(jj);
                if i==q(ii)
                    sim_cov(j+1:j+p(jj),j+1:j+p(jj)) = cluster_cov(kk)/2;
                else
                    sim_cov(j+1:j+p(jj),j+1:j+p(jj)) = cluster_cov(kk);
                end
            end
            sim_cov(logical(diag(ones(n,1))))=1;
            sim_data = mvnrnd(zeros(n,1),sim_cov,m);
            sim_data = zscore(sim_data);

            % display covariance matrix
            figure(fhs(1));
            ah=subplot(4,2,ct);
            imagesc(cov(sim_data));
            k = n-q(ii)*(p(jj)-1);
            title({sprintf('num features = %i, clusters = %i',n,q(ii));...
                sprintf('effective dimensions = %i',k)});
            axis('equal','tight');
            colormap(ah,egoalley);
            colorbar;
            caxis([-1 1]);
            drawnow;

            % plot connected components output
            out = decathlonConnCompDropN(sim_data);
            figure(fhs(2));
            ah=subplot(4,2,ct);
            imagesc(out);
            hold on;
            plot([k k],[0.5 size(out,1)+.5],'k--');
            ylabel('features dropped');
            xlabel('effective dimensionality');
            title('connected components heatmap');
            axis('equal','tight');
            colormap(ah,logjet_cmap);
            colorbar; caxis([0 1]);
            drawnow;
            
            % plot connected components output
            [~,bin_cts]=decathlonConnCompSweep(sim_data,200);
            figure(fhs(3));
            ah=subplot(4,2,ct);
            histogram(bin_cts,1:size(sim_data,2),'Normalization','probability');
            hold on;
            plot([k k]+0.5,[0 1],'k--');
            set(ah,'YLim',[0 1],'YTick',0:0.2:1);
            ylabel('probability');
            xlabel('effective dimensionality');
            drawnow;
        end
    end
end
