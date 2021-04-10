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


%% NEW REVISION ANALYSES

%%


D_z = load_decathlon_structs(pwd,'D_zscored_unfilled');
D_p = pair_decathlon_structs(D_z(1:2),'ImputeMode','none','CollapseMode','none');
D_cat = cat_decathlon_structs(D_p,'ImputeMode','none','CollapseMode','none');

f_cat = D_cat.fields;

f2 = D_z(2).fields;
[a,m,~] = parse_fieldnames(f2);
f2(~strcmpi(a,'Circadian')) = cellfun(@(a,m) sprintf('%s %s',a,m), ...
    a(~strcmpi(a,'Circadian')),m(~strcmpi(a,'Circadian')),'UniformOutput',false);
f2 = cellfun(@(s) ~any(strcmpi(f_cat,s)), f2);

f_cat = cat(1,f_cat,D_z(2).fields(f2));
[a,b,grp_idx] = group_apriori_fields(f_cat);
grp_idx = cat(1,grp_idx{:});

sim_data = NaN(size(D_cat.data,1),size(D_z(2).data,2));
sim_data(:,1:sum(~f2)) = D_cat.data;
sim_data(193:end,sum(~f2)+1:end) = D_z(2).data(:,f2);
sim_data = sim_data(:,grp_idx);
is_missing = isnan(sim_data);
sim_data(isnan(sim_data))=0;
cov_mat = cov(sim_data);
g_truth = mvnrnd(zeros(size(sim_data)),cov_mat);

%%

frac_missing = linspace(.05,.4,8);
frac_mask = sum(is_missing)./numel(is_missing);
nreps = 100;

% mean imputation
mean_impute_mse = zeros(nreps,numel(frac_missing));
r_impute_mean = zeros(nreps,numel(frac_missing));
r_true = corr(g_truth);
r_mask = upper_triangle_idx(size(r_true,1));
r_true = r_true(r_mask);
for i = 1:numel(frac_missing)
    frac = frac_missing(i);
    fprintf('frac. missing %i of %i\n',i,numel(frac_missing))
    for j=1:nreps
       
       p = randperm(size(is_missing,2));
       data = g_truth;
       mask = is_missing(:,p);
       [~,cut_off] = min(abs(frac-cumsum(frac_mask(p))));
       mask = mask & [true(1,cut_off) false(1,numel(frac_mask)-cut_off)];
       data(mask) = NaN;
       data = nanzscore(data);
       data(mask) = 0;
       
       % calculate error
       r = corr(data);
       r_impute_mean(j,i) = corr(r(r_mask),r_true);
       residuals = data(mask) - g_truth(mask);
       mean_impute_mse(j,i) = mean(residuals.^2);
    end
end


% als imputation
nreps = 100;
als_impute_mse = zeros(nreps,numel(frac_missing));
clust_als_impute_mse = zeros(nreps,numel(frac_missing));
non_clust_als_impute_mse = zeros(nreps,numel(frac_missing));
r_impute_als = zeros(nreps,numel(frac_missing));
als_data = cell(numel(frac_missing),1);
for i = 1:numel(frac_missing)
    frac = frac_missing(i);
    [~,cut_off] = min(abs(frac-cumsum(frac_mask)));
    parfor j=1:nreps
        fprintf('rep %i of %i\n',j,nreps)
        p = randperm(size(is_missing,2));
        D.data = g_truth;
        mask = is_missing;
        mask = mask & [true(1,cut_off) false(1,numel(frac_mask)-cut_off)];
        D.data(mask) = NaN;
        [D,~] = avg_als_impute(D,400);
        D.data = zscore(D.data);

        % calculate error
        r = corr(D.data);
        r_impute_als(j,i) = corr(r(r_mask),r_true);
        residuals = D.data(mask) - g_truth(mask);
        als_impute_mse(j,i) = mean(residuals.^2);
        residuals = D.data(mask & clust_feat) - g_truth(mask & clust_feat);
        clust_als_impute_mse(j,i) = mean(residuals.^2);
        residuals = D.data(mask & non_clust_feat) - g_truth(mask & non_clust_feat);
        non_clust_als_impute_mse(j,i) = mean(residuals.^2);
    end
    
end

%% Plot mean vs als imputation MSE relative to ground truth as a function of fraction missing


figure; hold on;
y = mean(mean_impute_mse);
ci95 = [prctile(mean_impute_mse,97.5);prctile(mean_impute_mse,2.5)];
vx = [frac_missing fliplr(frac_missing)];
vy = [ci95(2,:) fliplr(ci95(1,:))];
patch('XData',vx(:),'YData',vy(:),'FaceColor',[.75 .75 .75],'EdgeColor','none');
plot(frac_missing,y,'k--');
y = mean(clust_als_impute_mse,1);
ci95 = [prctile(clust_als_impute_mse,97.5,1);prctile(clust_als_impute_mse,2.5,1)];
vx = [frac_missing fliplr(frac_missing)];
vy = [ci95(2,:) fliplr(ci95(1,:))];
patch('XData',vx(:),'YData',vy(:),'FaceColor',[.75 .55 .5],'EdgeColor','none');
plot(frac_missing,y,'r');
set(gca,'XLim',frac_missing([1 end]),'YLim',[0 1]);
xlabel('fraction missing');
ylabel('MSE');

%% Plot ALS vs MEAN imputation as a function of frac missing for simulated
%  data generated on decathlon covariance matrix

frac_missing = linspace(.05,.4,8);
als_data = cell(numel(frac_missing),1);
frac_mask = sum(is_missing)./numel(is_missing);
parfor i = 1:numel(frac_missing)
    frac = frac_missing(i);
    [~,cut_off] = min(abs(frac-cumsum(frac_mask)));
    fprintf('fraction missing %i of %i\n',i,numel(frac_missing))
    D.data = g_truth;
    mask = is_missing;
    mask = mask & [true(1,cut_off) false(1,numel(frac_mask)-cut_off)];
    D.data(mask) = NaN;
    [D,~] = avg_als_impute(D,50);
    als_data{i} = zscore(D.data);
end
mean_data = cell(numel(frac_missing),1);
for i = 1:numel(frac_missing)
    frac = frac_missing(i);
    fprintf('frac. missing %i of %i\n',i,numel(frac_missing))
       
   p = randperm(size(is_missing,2));
   data = g_truth;
   mask = is_missing(:,p);
   [~,cut_off] = min(abs(frac-cumsum(frac_mask(p))));
   mask = mask & [true(1,cut_off) false(1,numel(frac_mask)-cut_off)];
   data(mask) = NaN;
   data = nanzscore(data);
   data(mask) = 0;
   mean_data{i} = data;
end


% cluster features
f1 = figure;
mean_cov = mean(abs(cov_mat));
high_cov = mean_cov > prctile(mean_cov,75);
low_cov = mean_cov < prctile(mean_cov,25);
ncols = ceil(sqrt(numel(frac_missing)));
nrows = ceil(numel(frac_missing)/ncols);
for i=1:numel(frac_missing)
    subplot(nrows,ncols,i); hold on;
    frac = frac_missing(i);
    [~,cut_off] = min(abs(frac-cumsum(frac_mask)));
    mask = is_missing;
    mask = mask & [true(1,cut_off) false(1,numel(frac_mask)-cut_off)];
    infilled = als_data{i}(mask & high_cov);
    gt = g_truth(mask & high_cov);
    pretty_scatter(gt,infilled,'r','MarkerSize',1.5); hold on;
    r1 = corr([gt infilled],'Type','pearson');

    % non cluster features
    infilled = als_data{i}(mask & low_cov);
    gt = g_truth(mask & low_cov);
    r2 = corr([gt infilled],'Type','pearson');
    pretty_scatter(gt,infilled,'k','MarkerSize',1.25);
    ylabel('imputed');
    xlabel('ground truth');
    axis('equal');
    legend({sprintf('high cov. r=%0.2f',r1(1,2));...
        sprintf('low-cov r=%0.2f',r2(1,2))});
    title(sprintf('fraction missing = %0.2f',frac));
    set(gca,'XLim',[-4 4],'YLim',[-4 4]);
end

% generate scree plots for ground truth, ALS, and mean imputed data
f2 = figure;
for i=1:numel(frac_missing)
    subplot(nrows,ncols,i); hold on;
    frac = frac_missing(i);
    plot_pca_bootstrap(mean_data{i},100,95,'noncummulative',[1 0 0]);
    plot_pca_bootstrap(als_data{i},100,95,'noncummulative',[0 1 1]);
    plot_pca_bootstrap(g_truth,100,95,'noncummulative',[0 0 1]);

    ah = gca;
    del = [1,2,6,7,9,10,12,14,15];
    lh = findall(ah,'Type','Line');
    arrayfun(@(h) delete(h), lh(del));
    ph = findall(ah,'Type','patch');
    del = [2,4];
    arrayfun(@(h) delete(h), ph(del));
    lh = findall(gca,'Type','Line');
    legend(lh([1,3,4,6]),{'observed mean';'shuffled mean';'als imputed';'mean imputed'});
    title(sprintf('fraction missing = %0.2f',frac));
    ah.YScale = 'log';
end

%% Plot simulated covariance matrix and sort by mean covariance

corr_mat = corr(sim_data);
mean_corr = mean(abs(corr_mat));
[~,p] = sort(mean_corr,'descend');
figure;
subplot(2,1,1);
imagesc(corr(sim_data));
colormap(nanticoke);
caxis([-1 1]);
axis('equal','tight','off');
colorbar;
set(gca,'TickLength',[0 0],'YTick',[],'XTick',[]);

subplot(2,1,2);

imagesc(corr_mat(p,p));
hold on;
colormap(nanticoke);
caxis([-1 1]);
axis('equal','tight','off');
colorbar;
set(gca,'TickLength',[0 0],'YTick',[],'XTick',[],'Clipping','off');
vx = [-1;-1;-.5;-.5;-1];
vy = floor(numel(mean_corr)*.25) .* [0;1;1;0;0];
vy = repmat(vy,1,2) + [0 ceil(numel(mean_corr)*.75)];
patch('XData',repmat(vx,1,2),'YData',vy,'FaceColor','r','EdgeColor','none');

