
save_dir = ['C:\Users\winsl0w\Documents\decathlon\decathlon_analysis\'...
    'matrices\decathlon_paper\rna_seq\lasso_models'];

% quantile normalize reads
reads = RNAseq.reads;
reads(:,sum(reads)<.5E6) = NaN;
reads = quantile_normalize(reads);
filt = find(~all(isnan(reads)));
x = reads(:,filt)';

% Perform modeling on FULL decathlon behavioral data

% load decathlon behavioral data
D = load_decathlon_structs(fdir,'D13_als');
opts = {'CollapseMode','PCA','CollapseFields','all'};
D_p = pair_decathlon_structs(D,opts{:});
d = NaN(size(D_p(1).data,1)+2,size(D_p(1).data,2));
mask = ~ismember(1:size(d,1),[99 114]);
d(mask,:) = D_p(1).data;
full_data = d(filt,:);
full_shuf_data = shuffle_columns(full_data);

nfields = size(d,2);
for i=1:nfields
    
    fprintf('iter %i of %i\n',i,nfields);
    
    % select behavioral metric to be modeled
    y_obs = full_data(:,i);
    y_shuf = full_shuf_data(:,i);
    
    % model observed data
    out_obs = decathlonLasso(x,y_obs);
    save(sprintf('%s\\lasso_coefs_obs_full_%i.mat',save_dir,i),'out_obs');
    clear out_obs
    
    % model shuffled data
    out_shuf = decathlonLasso(shuffle_columns(x),y_shuf);
    save(sprintf('%s\\lasso_coefs_shuf_full_%i.mat',save_dir,i),'out_shuf');
    clear out_shuf
end

% plot r-value distributions
n = size(full_data,1);
r_obs = arrayfun(@(s) nanmean(s.rs),obs_full_mdls);
r_shuf = arrayfun(@(s) nanmean(s.rs),shuf_full_mdls);
t_obs = abs(r_obs)./sqrt((1-abs(r_obs).^2)./(n-2));
p_obs = (1-tcdf(t_obs,n))*2;

autoPlotDist(r_obs,true(size(r_obs)));
autoPlotDist(r_shuf,true(size(r_obs)),gca);

%% Perform modeling on DISTILLED decathlon data

% rnaseq data
filt = find(~all(isnan(reads)));
x = reads(:,filt)';

% load decathlon behavioral data
D = load_decathlon_structs;
opts = {'CollapseMode','PCA','CollapseFields','all',...
    'Standardize',false,'ImputeMode','none'};
D_p = pair_decathlon_structs(D,opts{:});
d = NaN(size(D_p(1).data,1)+2,size(D_p(1).data,2));
mask = ~ismember(1:size(d,1),[99 114]);
d(mask,:) = D_p(1).data;
[distilled_pca_coefs,distilled_data] = pca(d(filt,:),'NumComponents',14);
distilled_data = zscore(distilled_data);
distilled_shuf_data = shuffle_columns(distilled_data);
nfields = size(distilled_data,2);
for i=1:nfields
    
    fprintf('iter %i of %i\n',i,nfields);
    
    % select behavioral metric to be modeled
    y_obs = distilled_data(:,i);
    y_shuf = distilled_shuf_data(:,i);
    
    % model observed data
    out_obs = decathlonLasso(x,y_obs);
    save(sprintf('%s\\lasso_coefs_obs_distilled_pcs_%i.mat',save_dir,i),'out_obs');
    clear out_obs
    
    % model shuffled data
    out_shuf = decathlonLasso(shuffle_columns(x),y_shuf);
    save(sprintf('%s\\lasso_coefs_shuf_distilled_pcs_%i.mat',save_dir,i),'out_shuf');
    clear out_shuf
end

% % plot r-value distributions
% n = size(distilled_data,1);
% r_obs = arrayfun(@(s) nanmean(s.rs),obs_distilled_mdls);
% r_shuf = arrayfun(@(s) nanmean(s.rs),shuf_distilled_mdls);
% t_obs = abs(r_obs)./sqrt((1-abs(r_obs).^2)./(n-2));
% p_obs = (1-tcdf(t_obs,n))*2;
% 
% figure;
% autoPlotDist(r_obs,true(size(r_obs)));
% autoPlotDist(r_shuf,true(size(r_obs)));


%%

obs_fpaths = recursiveSearch(save_dir,'keyword','obs_full');
shuf_fpaths = recursiveSearch(save_dir,'keyword','shuf_full');

obs_full_r = cell(numel(obs_fpaths),1);
shuf_full_r = cell(numel(obs_fpaths),1);
obs_n_coefs = cell(numel(obs_fpaths),1);
shuf_n_coefs = cell(numel(obs_fpaths),1);
for i=1:numel(obs_fpaths)
   fprintf('loading file %i of %i\n',i,numel(obs_fpaths));
   load(obs_fpaths{i},'out_obs');
   load(shuf_fpaths{i},'out_shuf');
   obs_full_r{i} = out_obs.rs;
   shuf_full_r{i} = out_shuf.rs;
    obs_n_coefs{i} = sum(out_obs.Bs~=0);
   shuf_n_coefs{i} = sum(out_shuf.Bs~=0);
   clear out_obs out_shuf
end

n=101;
obs_r = cat(2,obs_full_r{:});
t_obs = abs(obs_r(:))./sqrt((1-abs(obs_r(:)).^2)./(n-2));
p_obs = (1-tcdf(t_obs,n))*2;

shuf_r = cat(2,shuf_full_r{:});
t_shuf = abs(shuf_r(:))./sqrt((1-abs(shuf_r(:)).^2)./(n-2));
p_shuf = (1-tcdf(t_shuf,n))*2;

figure;
subplot(2,1,1);
imagesc(obs_r');
title('r-values - observed full data');
xlabel('resamplings');
ylabel('behavioral metrics');
colormap(egoalley); colorbar; caxis([-1 1]);

subplot(2,1,2);
imagesc(shuf_r');
title('r-values - shuffled full data');
xlabel('resamplings');
ylabel('behavioral metrics');
colormap(egoalley); colorbar; caxis([-1 1]);


obs_n_coefs = cat(1,obs_n_coefs{:});
shuf_n_coefs = cat(1,shuf_n_coefs{:});

subplot(1,3,1);
histogram(obs_r(:),linspace(-1,1,50),'Normalization','PDF'); hold on;
histogram(shuf_r(:),linspace(-1,1,50),'Normalization','PDF');
ylabel('PDF');
xlabel('r-value');
legend({'observed';'shuffled'});
title('D1 - GxB models (full)');

subplot(1,3,2);
histogram(p_obs,linspace(0,1,50),'Normalization','PDF'); hold on;
histogram(p_shuf,linspace(0,1,50),'Normalization','PDF');
ylabel('PDF');
xlabel('p-value');
legend({'observed';'shuffled'});
title('D1 - GxB models (full)');

subplot(1,3,3);
histogram(obs_n_coefs(:),linspace(0,100,100),'Normalization','PDF'); hold on;
histogram(shuf_n_coefs(:),linspace(0,100,100),'Normalization','PDF');
ylabel('PDF');
xlabel('lasso - num. coefficients');
legend({'observed';'shuffled'});
title('D1 - GxB models (full)');

%% distribution of number of coefficients

obs_fpaths = recursiveSearch(save_dir,'keyword','obs_distilled_pcs');
shuf_fpaths = recursiveSearch(save_dir,'keyword','shuf_distilled_pcs');

obs_distilled_r = cell(numel(obs_fpaths),1);
shuf_distilled_r = cell(numel(obs_fpaths),1);
obs_n_coefs = cell(numel(obs_fpaths),1);
shuf_n_coefs = cell(numel(obs_fpaths),1);
for i=1:numel(obs_fpaths)
   fprintf('loading file %i of %i\n',i,numel(obs_fpaths));
   load(obs_fpaths{i},'out_obs');
   load(shuf_fpaths{i},'out_shuf');
   obs_distilled_r{i} = out_obs.rs;
   shuf_distilled_r{i} = out_shuf.rs;
   obs_n_coefs{i} = sum(out_obs.Bs~=0);
   shuf_n_coefs{i} = sum(out_shuf.Bs~=0);
   clear out_obs out_shuf
end

n=101*0.1;
obs_r = cat(2,obs_distilled_r{:});
t_obs = abs(obs_r(:))./sqrt((1-abs(obs_r(:)).^2)./(n-2));
p_obs = (1-tcdf(t_obs,n))*2;

shuf_r = cat(2,shuf_distilled_r{:});
t_shuf = abs(shuf_r(:))./sqrt((1-abs(shuf_r(:)).^2)./(n-2));
p_shuf = (1-tcdf(t_shuf,n))*2;

obs_n_coefs = cat(1,obs_n_coefs{:});
shuf_n_coefs = cat(1,shuf_n_coefs{:});

figure;
subplot(2,1,1);
imagesc(obs_r');
title('r-values - observed distilled data');
xlabel('resamplings');
ylabel('behavioral metrics');
colormap(egoalley); colorbar; caxis([-1 1]);

subplot(2,1,2);
imagesc(shuf_r');
title('r-values - shuffled distilled data');
xlabel('resamplings');
ylabel('behavioral metrics');
colormap(egoalley); colorbar; caxis([-1 1]);

figure;
subplot(1,3,1);
histogram(obs_r(:),linspace(-1,1,50),'Normalization','PDF'); hold on;
histogram(shuf_r(:),linspace(-1,1,50),'Normalization','PDF');
ylabel('PDF');
xlabel('r-value');
legend({'observed';'shuffled'});
title('D1 - GxB models (distilled)');

subplot(1,3,2);
histogram(p_obs,linspace(0,1,50),'Normalization','PDF'); hold on;
histogram(p_shuf,linspace(0,1,50),'Normalization','PDF');
ylabel('PDF');
xlabel('p-value');
legend({'observed';'shuffled'});
title('D1 - GxB models (distilled)');

subplot(1,3,3);
histogram(obs_n_coefs(p_obs<0.05),linspace(0,100,100),'Normalization','PDF'); hold on;
histogram(shuf_n_coefs(p_shuf<0.05),linspace(0,100,100),'Normalization','PDF');
ylabel('PDF');
xlabel('lasso - num. coefficients');
legend({'observed';'shuffled'});
title('D1 - GxB models (distilled)');



%%

y = d(filt,:);

% quantile normalize reads
npcs = 8;
[~,x] = pca(reads(:,filt)','NumComponents',npcs);

nfields = size(d,2);
lin_rvals = cell(nfields,1);
lin_shuf_rvals = cell(nfields,1);
[~,y_obs] = pca(y,'NumComponents',14);
[~,y_shuf] = pca(shuffle_columns(y),'NumComponents',14);
nfields = size(y_obs,2);
for i=1:nfields
    fprintf('iter %i of %i\n',i,nfields);
    obs_mdl = decathlonLM(x,y_obs(:,i));
    lin_rvals{i} = obs_mdl.rs;
    shuf_mdl = decathlonLM(shuffle_columns(x),y_shuf(:,i));
    lin_shuf_rvals{i} = shuf_mdl.rs;
end

n=101;
obs_r = nanmean(cat(2,lin_rvals{:}));
t_obs = abs(obs_r(:))./sqrt((1-abs(obs_r(:)).^2)./(n-2));
p_obs = (1-tcdf(t_obs,n))*2;

shuf_r = nanmean(cat(2,lin_shuf_rvals{:}));
t_shuf = abs(shuf_r(:))./sqrt((1-abs(shuf_r(:)).^2)./(n-2));
p_shuf = (1-tcdf(t_shuf,n))*2;

figure;
subplot(1,2,1);
histogram(obs_r(:),linspace(-1,1,50),'Normalization','PDF'); hold on;
histogram(shuf_r(:),linspace(-1,1,50),'Normalization','PDF');
ylabel('PDF');
xlabel('r-value');
legend({'observed';'shuffled'});
title('D1 - GxB models (distilled)');

subplot(1,2,2);
histogram(p_obs,linspace(0,1,50),'Normalization','PDF'); hold on;
histogram(p_shuf,linspace(0,1,50),'Normalization','PDF');
ylabel('PDF');
xlabel('p-value');
legend({'observed';'shuffled'});
title('D1 - GxB models (distilled)');

figure;
subplot(1,2,1);
bar(-log(p_obs));
hold on; plot([0 size(x,2)+1],-log([0.05 0.05]),'k--');
ylabel('-log(p)');
xlabel('distilled matrix PCs');
title('observed data');
set(gca,'YLim',[0 7]);

subplot(1,2,2);
bar(-log(p_shuf));
hold on; plot([0 size(x,2)+1],-log([0.05 0.05]),'k--');
ylabel('-log(p)');
xlabel('distilled matrix PCs');
set(gca,'YLim',[0 7]);

title('shuffled data');

%%

% quantile normalize reads
reads = quantile_normalize(RNAseq.reads);
filt = find(~all(isnan(reads)));
npcs = 14;
%x = zscore(reads(:,filt)');
x = reads(:,filt)';
gene_variance = var(x);
gene_mean = mean(x);
[pca_coefs,x] = pca(x,'NumComponents',8);
x=zscore(x);

% load decathlon behavioral data
D = load_decathlon_structs;
opts = {'CollapseMode','PCA','CollapseFields','all',...
    'Standardize',false,'ImputeMode','none'};
D_p = pair_decathlon_structs(D,opts{:});
d = NaN(size(D_p(1).data,1)+2,size(D_p(1).data,2));
mask = ~ismember(1:size(d,1),[99 114]);
d(mask,:) = D_p(1).data;
distilled_data = d(filt,:);
rng('shuffle');
distilled_shuf_data = shuffle_columns(distilled_data);

[beh_pc_loadings,y_obs] = pca(distilled_data,'NumComponents',npcs);
y_shuf = shuffle_columns(y_obs);

nfields = size(y_obs,2);
lin_obs_rvals = cell(nfields,1);
lin_shuf_rvals = cell(nfields,1);
lin_obs_gene_coefs = cell(nfields,1);
lin_obs_gene_coefs(:) = {pca_coefs};
lin_shuf_gene_coefs = cell(nfields,1);
lin_shuf_gene_coefs(:) = {pca_coefs};

for i=1:nfields
    fprintf('iter %i of %i\n',i,nfields);
    obs_mdl = decathlonLM(x,y_obs(:,i));
    lin_obs_rvals{i} = obs_mdl.rs;
    lin_obs_gene_coefs{i} = sum(lin_obs_gene_coefs{i}.*...
        nanmedian(cat(2,obs_mdl.Bs{:}),2)',2);
    shuf_mdl = decathlonLM(x,y_shuf(:,i));
    lin_shuf_rvals{i} = shuf_mdl.rs;
    lin_shuf_gene_coefs{i} = sum(lin_shuf_gene_coefs{i}.*...
        nanmedian(cat(2,shuf_mdl.Bs{:}),2)',2);
end

n=size(x,1)-1;
obs_r = nanmean(cat(2,lin_obs_rvals{:}),1);
t_obs = abs(obs_r(:))./sqrt((1-abs(obs_r(:)).^2)./(n-2));
p_obs = (1-tcdf(t_obs,n))*2;

shuf_r = nanmean(cat(2,lin_shuf_rvals{:}),1);
t_shuf = abs(shuf_r(:))./sqrt((1-abs(shuf_r(:)).^2)./(n-2));
p_shuf = (1-tcdf(t_shuf,n))*2;

figure;
subplot(1,2,1);
bar(-log10(p_obs));
hold on; plot([0 size(y_obs,2)+1],-log10([0.05 0.05]),'k--');
ylabel('-log(p)');
xlabel('distilled matrix PCs');
title('observed data');
set(gca,'YLim',[0 3]);

subplot(1,2,2);
bar(-log10(p_shuf));
hold on; plot([0 size(y_obs,2)+1],-log10([0.05 0.05]),'k--');
ylabel('-log(p)');
xlabel('distilled matrix PCs');
set(gca,'YLim',[0 3]);
title('shuffled data');

% plot top loadings for significant models
[~,sig_models] = sort(p_obs);
sig_models = sig_models(1:3);
lin_candidate_ids = cell(numel(sig_models),1);
lin_shuf_ids = cell(numel(sig_models),1);
lin_candidate_variance = cell(numel(sig_models),1);
lin_shuf_variance = cell(numel(sig_models),1);


%%
out_path = ['C:\Users\winsl0w\Documents\decathlon\decathlon_analysis\'...
    'matrices\decathlon_paper\rna_seq\gene_ids\'];
n_max = numel(sig_models);
for i=1:n_max
   
    idx = sig_models(i);
    coefs =  lin_obs_gene_coefs{idx};
    coefs= abs(coefs);
    coef_cdf = cumsum(coefs./sum(coefs));
    lb = prctile(coefs,95);
    lin_candidate_ids{i} = RNAseq.geneID(coefs > lb);
    lin_candidate_variance{i} = gene_variance(coefs > lb);

    fid = fopen(sprintf('%slinearmodel_geneids_distilledpc_%i.txt',out_path,idx),'W+');
    out_txt = cellfun(@(s) sprintf('%s\n',s), lin_candidate_ids{i}, 'UniformOutput', false);
    fprintf(fid,'%s',cat(2,out_txt{:}));
    fclose(fid);

    coefs =  lin_obs_gene_coefs{idx};
    [coefs,p] = sort(coefs);
    ids = RNAseq.geneID(p);
    ids = [ids(1:10);ids(end-9:end)];
    coefs = [coefs(1:10);coefs(end-9:end)];

    % make plot
    subplot(1,n_max,i);
    barh(coefs);
    title(sprintf('distilled behavior PC %i',idx));
    set(gca,'YTick',1:numel(ids),'YTickLabel',ids,'TickLength',[0 0],...
        'FontSize',6);

    coefs =  lin_shuf_gene_coefs{idx};
    lb = prctile(coefs,2.5);
    ub = prctile(coefs,97.5);
    lin_shuf_ids{i} = RNAseq.geneID(coefs < lb | coefs > ub);
    lin_shuf_variance{i} = gene_variance(coefs < lb | coefs > ub);
end


candidates = false(n_max,size(reads,1));
for i=1:n_max
    
    idx = sig_models(i);
    coefs =  lin_obs_gene_coefs{idx};
    lb = prctile(coefs,2.5);
    ub = prctile(coefs,97.5);
    candidates(i,:) = coefs < lb | coefs > ub;
end

figure; hold on;
colors = {[0 0 0];[0 0 0.9];[0.1 0.9 0.1];[0.9 0 0];[1 0.5 0];[1 0.5 1]};
candidate_cts = sum(candidates);
unique_cts = unique(candidate_cts);
for i=1:numel(unique_cts)
    
    f = candidate_cts == unique_cts(i);
    vx = gene_mean(f);
    vy = gene_variance(f);
    ph = pretty_scatter(log(vx),log(vy),colors{i},'MarkerSize',4);
    xlabel('log mean expression');
    ylabel('log expression variance');
    title(sprintf('linear model %i',i));
end

n_max = numel(sig_models);
for i=1:n_max
   
    idx = sig_models(i);
    tmp_can = candidates(i,:) & candidate_cts < 2;
    ids = RNAseq.geneID(tmp_can);

    fid = fopen(sprintf('%sunique_ids_distpc_%i.txt',out_path,idx),'W+');
    out_txt = cellfun(@(s) sprintf('%s\n',s), ids, 'UniformOutput', false);
    fprintf(fid,'%s',cat(2,out_txt{:}));
    fclose(fid);
end

figure;
for i=1:n_max
   % make plot
    subplot(1,n_max,i);
    histogram(log(lin_shuf_variance{i}));
end

% plot loadings for significant behavioral PCs
figure;
for i=1:n_max
    idx = sig_models(i);
    [coefs,p] = sort(beh_pc_loadings(:,idx));
    
    subplot(1,n_max,i);
    barh(coefs);
    set(gca,'YTick',1:numel(coefs),'YTickLabel',D_p(1).fields(p),...
        'TickLength',[0 0]);
    title(sprintf('distilled PC %i',idx));
end


%% construct gene lists for the venn diagram

can_pairs = cell(7,1);
can_pairs(1:3) = num2cell(1:3);
can_pairs(4:6) = num2cell(unique_idx_pairs(3,1),2);
can_pairs(7) = {1:3};
can_mask = cell(7,1);
for i=1:numel(can_pairs)
    include = can_pairs{i};
    exclude = find(~ismember(1:3,include));
    mask = [candidates(include,:);~candidates(exclude,:)];
    can_mask{i} = sum(mask,1)==3;
end

for i=1:numel(can_mask)
    pc_str = sprintf(repmat('-%i',1,numel(can_pairs{i})),sig_models(can_pairs{i}));
    fid = fopen(sprintf('%svenn_geneids_pcs%s_(%i).txt',out_path,pc_str,sum(can_mask{i})),'W+');
    ids = RNAseq.geneID(can_mask{i});
    out_txt = cellfun(@(s) sprintf('%s\n',s), ids, 'UniformOutput', false);
    fprintf(fid,'%s',cat(2,out_txt{:}));
    fclose(fid);
end

figure;
n_genes = cellfun(@sum,can_mask);
[~,vh] = venn(n_genes(1:3),n_genes(4:end));

% label zone centers
for i=1:numel(sig_models)
    cen = vh.Position(i,:);
    text(cen(1),cen(2),sprintf('PC%i',sig_models(i)),'HorizontalAlignment','center');
end
for i=1:numel(can_mask)
    cen = vh.ZoneCentroid(i,:);
    text(cen(1),cen(2),sprintf('n=%i',n_genes(i)),'HorizontalAlignment','center');
end
axis('equal','off');

%% plot heatmaps of canoncorr p-values as a function of no. pcs of gene expression

k_to_try = 8:30;
ah = subplot_array(numel(k_to_try)*2);
n = size(d,2);

for j=1:numel(k_to_try)
    k=k_to_try(j);
    x = reads(:,filt)';
    [~,x] = pca(x,'NumComponents',k);
%     [~,y_obs] = pca(distilled_data,'NumComponents',n);
%     [~,y_shuf] = pca(distilled_shuf_data,'NumComponents',n);
    y_obs = distilled_data;
    y_shuf = distilled_shuf_data;
    canon_val_p = zeros(size(x,2)-size(y_obs,2),10);
    shuf_val_p = zeros(size(x,2)-size(y_obs,2),10);

    ct=1;
    for i=k:size(y_obs,2)
        [~,~,rnd_canon_r,~,~,rnd_canon_stats] = canoncorr(x,y_obs(:,1:i));
        canon_val_p(ct,1:numel(rnd_canon_stats.p)) = rnd_canon_stats.p;
        [~,~,rnd_canon_r,~,~,rnd_canon_stats] = canoncorr(x,y_shuf(:,1:i));
        shuf_val_p(ct,1:numel(rnd_canon_stats.p)) = rnd_canon_stats.p;
        ct = ct+1;
    end

    imagesc(canon_val_p,'Parent',ah((j-1)*2+1));
    ylabel('No. predictors in X');
    xlabel('canon vars');
    cb = colorbar(ah((j-1)*2+1));
    cb.Label.String = 'p-value';
    title(ah((j-1)*2+1),'observed');
    
    imagesc(shuf_val_p,'Parent',ah((j-1)*2+2));
    ylabel('No. predictors in X');
    xlabel('canon vars');
    cb = colorbar(ah((j-1)*2+2));
    cb.Label.String = 'p-value';
    title(ah((j-1)*2+2),'shuffled');
end

%% perform canonical correlations analysis on the gene pcs x distilled metrics


x = zscore(reads(:,filt)');
[gene_pca_coefs,x,gene_lat,~,gene_pc_vexp] = pca(x,'NumComponents',8);
gene_lat = gene_lat(1:size(x,2));
y_obs = nanzscore(distilled_data);
y_shuf = nanzscore(distilled_shuf_data);


[obs_a,obs_b,obs_canon_r,obs_u,obs_v,obs_canon_stats] = canoncorr(x,y_obs);
[~,~,shuf_canon_r,~,~,shuf_canon_stats] = canoncorr(x,y_shuf);

figure; hold on;
plot(-log10(obs_canon_stats.p),'r-');
plot(-log10(shuf_canon_stats.p),'k-');
plot([1 numel(shuf_canon_r)],-log10([0.05 0.05]),'k--');

% unpack behaviors in significant canonical variables
max_n = 3;
alpha = 5;
out_path = ['C:\Users\winsl0w\Documents\decathlon\decathlon_analysis\'...
    'matrices\decathlon_paper\rna_seq\gene_ids\'];
candidates = cell(max_n,1);
canon_candidate_coefs = cell(max_n,1);
canon_candidate_ids = cell(max_n,1);
canon_candidate_variance = cell(max_n,1);

for i=1:max_n
    subplot(1,max_n,i)
    [canon_beh_coefs,p] = sort(obs_b(:,i));
    barh(canon_beh_coefs);
    set(gca,'YTick',1:numel(p),'YTickLabel',D_p(1).fields(p),'TickLength',[0 0],...
        'FontSize',6);
    
    [canon_gene_coefs,p] = sort(obs_a(:,i));
    comb_gene_coefs = sum(gene_pca_coefs.*canon_gene_coefs',2);
    lower_bound = prctile(comb_gene_coefs,alpha/2);
    upper_bound = prctile(comb_gene_coefs,100-alpha/2);
    candidates{i} = comb_gene_coefs < lower_bound | comb_gene_coefs > upper_bound;
    canon_candidate_coefs{i}  = comb_gene_coefs(candidates{i});
    canon_candidate_ids{i} = RNAseq.geneID(candidates{i});
    canon_candidate_variance{i} = gene_variance(candidates{i});
    
    fid = fopen(sprintf('%scanoncorr_geneids_%i.txt',out_path,i),'W+');
    out_txt = cellfun(@(s) sprintf('%s\n',s), canon_candidate_ids{i}, 'UniformOutput', false);
    fprintf(fid,'%s',cat(2,out_txt{:}));
    fclose(fid);
    
    % get full behavioral matrix coefficients
    grps = regexp(D_p(1).fields,' ','split');
    grps = cellfun(@(s) s{1}, grps, 'UniformOutput', false);
    [~,grp_bounds,grps] = unique(grps,'stable');

    for j=1:numel(grp_bounds)
        if j==numel(grp_bounds)
            idx = grp_bounds(j):numel(grp_bounds);
        else
            idx = grp_bounds(j):grp_bounds(j+1);
        end
        
    end
end

figure;
for i=1:max_n
    subplot(1,max_n,i);
    [coefs,p] = sort(canon_candidate_coefs{i});
    ids = canon_candidate_ids{i}(p);
    ids = [ids(1:10);ids(end-9:end)];
    coefs = [coefs(1:10);coefs(end-9:end)];
    barh(coefs);
    title(sprintf('canon. var. %i',i));
    set(gca,'YTick',1:numel(ids),'YTickLabel',ids,'TickLength',[0 0],...
        'FontSize',6);
end

figure;
for i=1:n_max
   % make plot
    subplot(1,n_max,i);
    histogram(log(canon_candidate_variance{i}));
end

%% compute overlap in candidate genes identified by PCA and canoncorr

num_lin = numel(lin_candidate_ids);
num_canon = numel(canon_candidate_ids);
gene_id_overlap = NaN(num_lin,num_canon);
for i=1:num_lin
    for j=1:num_canon
        gene_id_overlap(i,j) = sum(ismember(lin_candidate_ids{i},...
            canon_candidate_ids{j}))/numel(lin_candidate_ids{i});
    end
end

figure;
subplot(2,1,1);
bh = bar(gene_id_overlap);
set(gca,'XTick',1:num_lin,'YLim',[0 1],'XTickLabel',...
    arrayfun(@(i) sprintf('linear model %i',i), 1:num_lin, 'UniformOutput', false));
legend(bh,arrayfun(@(i) sprintf('canon var %i',i), 1:num_lin, 'UniformOutput', false));
title('(observed) fraction overlap in significant gene lists');


num_lin = numel(lin_shuf_ids);
num_canon = numel(canon_candidate_ids);
gene_id_overlap = NaN(num_lin,num_canon);
for i=1:num_lin
    for j=1:num_canon
        gene_id_overlap(i,j) = sum(ismember(lin_shuf_ids{i},...
            canon_candidate_ids{j}))/numel(lin_shuf_ids{i});
    end
end

subplot(2,1,2);
bh = bar(gene_id_overlap);
set(gca,'XTick',1:num_lin,'YLim',[0 1],'XTickLabel',...
    arrayfun(@(i) sprintf('linear model %i',i), 1:num_lin, 'UniformOutput', false));
legend(bh,arrayfun(@(i) sprintf('canon var %i',i), 1:num_lin, 'UniformOutput', false));
title('(shuffled) fraction overlap in significant gene lists');

