clearvars -except RNAseq
load(['D:\decathlon_preprint_code_data_figures\'...
    'decathlon_paper_data\decathlon_rnaseq_data_2.mat']);
RNAseq = RNAseq(1);

% load decathlon data
D = load_decathlon_structs;

% generate distilled matrix and restrict to BerlinK-iso decathlon-2
opts = {'CollapseMode','PCA','CollapseFields','all',...
    'Standardize',true,'ImputeMode','none'};
D = pair_decathlon_structs(D,opts{:});
D = standardize_by_field(D);
D = D(1);
d = D.data(191:end,:);

% quantile normalize the reads
tmp_reads = RNAseq.raw_reads;
low_read_filt = nansum(tmp_reads) < 1E6;
fly_ids = RNAseq.fly_ids(~low_read_filt)-192;
tmp_reads = tmp_reads(:,~low_read_filt);
tmp_reads = quantile_normalize(tmp_reads);

norm_reads = NaN(size(tmp_reads,1),size(d,1));
norm_reads(:,fly_ids) = tmp_reads;

is_sequenced = any(~isnan(norm_reads));
norm_reads = norm_reads(:,is_sequenced);
d = d(is_sequenced,:);

%% bootstrap linear models

nreps = 1000;
[~,x] = pca(norm_reads','NumComponents',8);
[~,dist_pcs] = pca(d,'NumComponents',14);

coeffs_obs = cell(size(x,2),1);
r_obs = NaN(size(x,2),1);
p_obs = NaN(size(x,2),1);
p_null = NaN(nreps,1);
for i=1:size(dist_pcs,2)
    
    fprintf('Bootstrapping model %i of %i\n',i,size(d,2));
    
    % bootstrap observed data
    %y = d(:,i);
    y = dist_pcs(:,i);
    [coeffs_obs{i},r_obs(i),p_obs(i)] = bootstrap_linmodel(x,y);
end

% % generate null distribution of p-vals
% for i=1:nreps
%     if ~mod(i,100)
%         fprintf('%i\n',i);
%     end
%     x_shuf = shuffle_columns(x);
% %     [~,~,p_null(i)] = ...
% %         bootstrap_linmodel(x_shuf,d(:,randi(size(d,2),1)));
%     [~,~,p_null(i)] = ...
%         bootstrap_linmodel(x_shuf,dist_pcs(:,randi(size(dist_pcs,2),1)));
% end

% plot results
figure;
p_null = cellfun(@(p) p(:), p_null, 'UniformOutput', false);
p_null = cat(1,p_null{:});
pn = cell(numel(p_obs)+1,1);
pn(:) = {NaN(10,1)};
pn(end) = {-log10(p_null)};
%labels = [behavior.distilled_fields; {'Null'}];
labels = [D(1).fields; {'Null'}];
%labels = [arrayfun(@(i) sprintf('PC%i',i),1:numel(p_obs,2),'UniformOutput',false) {'Null'}]';
violinPlot(pn,'Labels',labels);
hold on;
bar([-log10(p_obs); NaN],0.5,'FaceColor',[.55 .55 .85],'EdgeColor','none');
set(gca,'YLim',[0 ceil(max(-log10(p_obs)))],'TIckLength',[0 0]);
plot([0 numel(p_obs)+2],[-log10(0.05) -log10(0.05)],'k--');
ylabel('-log(p)');


%%

x = norm_reads';
nreps = 1000;

coeffs_obs = cell(size(x,2),1);
p_obs = NaN(size(x,2),1);
r_obs = p_obs;
p_null = p_obs;
for i=1:size(d,2)
    
    fprintf('Bootstrapping model %i of %i\n',i,size(d,2));
    
    % bootstrap observed data
    y = d(:,i);
    %y = dist_pcs(:,i);
    out = decathlonLasso(x,y);
    p_obs(i) = out.ps;
    r_obs(i) = out.rs;
    coeffs_obs{i} = out.Bs;
end

% generate null distribution of p-vals
for i=1:nreps
    if ~mod(i,100)
        fprintf('%i\n',i);
    end
    x_shuf = shuffle_columns(x);
    out = decathlonLasso(x_shuf,d(:,randi(size(d,2),1)));
    p_null(i) = out.ps;
end

% plot results
figure;
pn = cell(numel(p_obs)+1,1);
pn(:) = {NaN(10,1)};
pn(end) = {p_null};
labels = [D(1).fields; {'Null'}];
violinPlot(pn,'Labels',labels);
hold on;
bar([-log10(p_obs); NaN],0.5,'FaceColor',[.55 .55 .85],'EdgeColor','none');
set(gca,'YLim',[0 ceil(max(-log10(p_obs)))],'TIckLength',[0 0]);
plot([0 numel(p_obs)+2],[-log10(0.05) -log10(0.05)],'k--');
ylabel('-log(p)');


%% perform canonical correlations analysis on the gene pcs x distilled metrics

D_p = D;
%x = zscore(norm_reads');
x=norm_reads';
gene_variance = var(norm_reads');
[gene_pca_coefs,x,gene_lat,~,gene_pc_vexp] = pca(x,'NumComponents',40);
gene_lat = gene_lat(1:size(x,2));
y_obs = nanzscore(d);
y_shuf = nanzscore(shuffle_columns(d));


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
    canon_candidate_ids{i} = RNAseq.geneIDs(candidates{i});
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
for i=1:max_n
   % make plot
    subplot(1,max_n,i);
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