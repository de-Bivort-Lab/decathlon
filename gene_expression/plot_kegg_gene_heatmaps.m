%% read in enrichment bootstrapped pvalues output by:
% 'shuffle_kegg_individual_gene_breakdown_(mean_log_p).r'

fdir = uigetdir(pwd,'Select RNAseq bootstrap output directory');

%%

shuf_p_paths = ...
    recursiveSearch(fdir,'ext','.csv','keyword','_shuffled_all_pval_mat');
unshuf_p_paths = ...
    recursiveSearch(fdir,'ext','.csv','keyword','_unshuffled_all_pval_mat');
shuf_n_paths = ...
    recursiveSearch(fdir,'ext','.csv','keyword','_shuffled_all_num_metrics_mat');
unshuf_n_paths = ...
    recursiveSearch(fdir,'ext','.csv','keyword','_unshuffled_all_num_metrics_mat');

% get metric matrices
shuf_metric_paths = ...
    recursiveSearch(fdir,'ext','.csv','keyword','_shuffled_all_cat_x_metric_mat');
unshuf_metric_paths = ...
    recursiveSearch(fdir,'ext','.csv','keyword','_unshuffled_all_cat_x_metric_mat');

d_p_shuf = repmat(struct('data',[],'rows',[],'cols',[]),numel(shuf_p_paths),1);
d_p_unshuf = d_p_shuf;
d_n_shuf = d_p_shuf;
d_n_unshuf = d_p_shuf;
d_m_shuf = d_p_shuf;
d_m_unshuf = d_p_shuf;

% abbreviate row labels
key_find = {'metabolism';'biosynthesis';'degradation';'alanine';'cystine';'serine';'proline';...
    'methionine';'tryptophan';'tyrosine';'lysine';'leucine';'glycine';'aspartate';...
    'glutamine';'valine';'threonine';'phenylalanine';'pathway';'fly';'-';...
    'autophagy animal';'cytochrome p450';'cysteine'};
key_replace = {'metab';'biosynth';'degrad';'ala';'cys';'ser';'pro';'met';'trp';'tyr';...
    'lys';'leu';'gly';'asp';'glu';'val';'thr';'phe';'';'';'';'autophagy';'CYPs';'cys'};


for i=1:numel(shuf_p_paths)
    
    % read in bootsrapped unshuffled pvals and num metrics 
    d_p_unshuf(i) = read_csv_mat(unshuf_p_paths{i});
    d_n_unshuf(i) = read_csv_mat(unshuf_n_paths{i});
    d_m_unshuf(i) = read_csv_mat(unshuf_metric_paths{i});


    % sort pvals, num metrics and enrichment cats by pvalue
    d_p_unshuf(i).data = nanmean(-log10(d_p_unshuf(i).data),2);
    [~,perm] = sort(d_p_unshuf(i).data);
    d_p_unshuf(i).rows = d_p_unshuf(i).rows(perm);
    d_p_unshuf(i).data = d_p_unshuf(i).data(perm);
    d_n_unshuf(i).data = nanmean(d_n_unshuf(i).data(perm,:),2);
    d_m_unshuf(i).rows = d_m_unshuf(i).rows(perm);
    d_m_unshuf(i).data = d_m_unshuf(i).data(perm,:);


    % read in bootsrapped shuffled pvals and num metrics 
    d_p_shuf(i) = read_csv_mat(shuf_p_paths{i});
    d_n_shuf(i) = read_csv_mat(shuf_n_paths{i});
    d_m_shuf(i) = read_csv_mat(shuf_metric_paths{i});
    
    % sort shuffled vals
    perm = cellfun(@(s) find(strcmp(d_p_shuf(i).rows,s)), d_p_unshuf(i).rows);
    d_p_shuf(i).rows = d_p_shuf(i).rows(perm);
    d_p_shuf(i).data = nanmean(-log10(d_p_shuf(i).data(perm,:)),2);
    d_n_shuf(i).data = nanmean(d_n_shuf(i).data(perm,:),2);
    d_m_shuf(i).rows = d_m_shuf(i).rows(perm);
    d_m_shuf(i).data = d_m_shuf(i).data(perm,:);
end

% change row labels
for i=1:numel(d_p_shuf)
    d_p_unshuf(i).rows = cellfun(@(s) lower(s), d_p_unshuf(i).rows, 'UniformOutput', false);
end

for k=1:numel(d_p_shuf)
    for i=1:numel(d_p_unshuf(k).rows)
        tmp_cat = d_p_unshuf(k).rows{i};
        for j=1:numel(key_find)
            start_idx = strfind(tmp_cat,key_find{j});
            if start_idx
                tmp_cat(start_idx:start_idx+numel(key_find{j})-1) = '';
                tmp_cat = cat(2,tmp_cat(1:start_idx-1),key_replace{j},tmp_cat(start_idx:end));
            end
        end
        d_p_unshuf(k).rows{i} = strtrim(tmp_cat);
    end
    d_m_unshuf(k).rows = d_p_unshuf(k).rows;
end


%% generate p-value bar plots for shuffled and unshuffled data

% enrichment p-value unshuffled/shuffled barplots
figure;
for i=1:numel(d_p_unshuf)
    subplot(numel(d_p_unshuf),1,i);  hold on;
    h1 = barh(d_p_unshuf(i).data,'FaceColor',[.6 .6 .6], 'EdgeColor', 'none');
    h2 = barh(d_p_shuf(i).data,'FaceColor',[0 0 0], 'EdgeColor', 'none');   
    set(gca,'YTick',1:numel(d_p_shuf(i).data),...
        'YTickLabel',d_p_unshuf(i).rows,'TickLength',[0 0],'XLim',[0 15]);
    legend([h1;h2],{'unshuffled','shuffled'},'Location','SouthEast');
    title(sprintf('decathlon %i',i));
    xlabel('-log_{10}[p-value]');
    ylabel('enrichment category');
    ax = gca;
    ax.YAxis.FontSize = 6;
    ax.YAxis.Label.FontSize = 10;
end
    

% enrichment p-value x num metrics scatter plot
figure;
for i=1:numel(d_p_unshuf)
    
    % scatter unshuffled data
    subplot(numel(d_p_unshuf),1,i);  hold on;
    lh1 = pretty_scatter(d_n_unshuf(i).data,d_p_unshuf(i).data,[.6 .6 .6],'MarkerSize',3);
    title(sprintf('decathlon %i',i));
    xlabel('num metrics');
    ylabel('-log_{10}[p]');
    
    % scatter shuffled data
    lh2 = pretty_scatter(d_n_shuf(i).data,d_p_shuf(i).data,[0 0 0],'MarkerSize',3);
    title(sprintf('decathlon %i',i));
    xlabel('num metrics');
    ylabel('-log_{10}[p]');
    set(gca,'XLim',[0 12],'YLim',[0 15]);
    legend([lh1;lh2],{'unshuffled','shuffled'},'Location','NorthWest');
    
    % apply text labels to significant nodes
    p_filt = find(d_p_unshuf(i).data>2);
    x = d_n_unshuf(i).data;
    y = d_p_unshuf(i).data;
    for j=1:numel(p_filt)
        text(x(p_filt(j))-.2,y(p_filt(j)),d_p_unshuf(i).rows(p_filt(j)),...
            'HorizontalAlignment','right','FontSize',6);
    end
end

% pair dot plot
[p_a,p_b] = ismember(d_p_unshuf(1).rows,d_p_unshuf(2).rows);
a = d_p_unshuf(1).data(p_a);
b = d_p_unshuf(2).data(p_b(p_b>0));
d2_lab = d_p_unshuf(2).rows(p_b(p_b>0));

figure;
subplot(1,2,1);
plot_pair_dot(a,b,'k',[]);
set(gca,'XTick',[1 2],'XTickLabel',{'D1';'D2'},'YLim',[0 15]);
title('unshuffled');
ylabel('-log_{10}[p-value]');
hold on;
for i=1:numel(d2_lab)
   text(2.1,b(i),d2_lab{i},'FontSize',6); 
end

[p_a,p_b] = ismember(d_p_shuf(1).rows,d_p_shuf(2).rows);
a = d_p_shuf(1).data(p_a);
b = d_p_shuf(2).data(p_b(p_b>0));


subplot(1,2,2);
plot_pair_dot(a,b,'k',[]);
set(gca,'XTick',[1 2],'XTickLabel',{'D1';'D2'});
title('shuffled');
ylabel('-log_{10}[p-value]');
set(gca,'YLim',[0 15]);

%% map metric names to row labels

% load behavioral data
fdir = uigetdir(pwd,'Select the decathlon_paper_data directory');
D_b = load_decathlon_structs(fdir,'D_als_filled');
D_b = D_b(1:2);
D_b = pair_decathlon_structs(D_b,'CollapseMode','PCA','CollapseFields','none');
fields = D_b(1).fields;
[~, ~, grp_idx] = group_apriori_fields(D_b(1));
grp_idx = cat(1,grp_idx{:});
for i=1:numel(d_m_unshuf)
    d_m_unshuf(i).cols = fields(grp_idx);
end

%% create a common sort/mapping for D1/D2 behavior x pathway heatmaps

% unique rows
u_rows = unique(cat(1,d_m_unshuf.rows));
u_cols = unique(cat(1,d_m_unshuf.cols));
d_sorted.data = zeros(numel(u_rows),numel(u_cols),numel(d_m_unshuf));
d_sorted.rows = cell(numel(u_rows),1);
d_sorted.cols = cell(numel(u_cols),1);
rows_shared = arrayfun(@(d) ...
    cellfun(@(s) any(strcmp(d.rows,s)), u_rows),...
    d_m_unshuf, 'UniformOutput', false);
rows_shared = u_rows(all(cat(2,rows_shared{:})'));

% sort rows common to all based on first dataset
idx_to_sort = find(cellfun(@(s) any(strcmp(s,rows_shared)), d_m_unshuf(1).rows));
[p_row,p_col] = get_cluster_perm(d_m_unshuf(1).data(idx_to_sort,:),'complete','cityblock');
sorted_rows = flip(d_m_unshuf(1).rows(idx_to_sort(p_row)));
sorted_cols = d_m_unshuf(1).cols(p_col);

for i=1:numel(d_m_unshuf)
    d = d_m_unshuf(i);
    unsorted_row_idx = find(~ismember(d.rows,sorted_rows));
    p_row = get_cluster_perm(d.data(unsorted_row_idx,:),'complete','cityblock');
    sorted_rows = cat(1,sorted_rows,flip(d.rows(unsorted_row_idx(p_row))));
    
    p_row = cellfun(@(s) find(strcmp(d.rows,s)), sorted_rows,'UniformOutput',false);
    sorted_row_mask = ~cellfun(@isempty,p_row);
    d_sorted.data(sorted_row_mask,:,i) = d.data(cat(1,p_row{:}),p_col);
end
d_sorted.rows = sorted_rows;
d_sorted.cols = sorted_cols;

figure;
for i=1:numel(d_m_unshuf)
    subplot(1,numel(d_m_unshuf),i);
    imagesc(fliplr(d_sorted.data(:,:,i)'));
    colormap(flip(bone));
    caxis([0 1]);
    cb = colorbar;
    cb.Label.String = 'Bootstrap probability';
    title(sprintf('decathlon batch #%i',i));
    set(gca,'XTick',1:numel(d_sorted.rows),'XTickLabels',...
        flip(d_sorted.rows),'XTickLabelRotation',90,...
        'YTick',1:numel(d_sorted.cols),'YTickLabels',...
        pretty_labels(d_sorted.cols),'FontSize',8,'TickLength',[0 0]);
end

%%

% plot matrices
figure;
for i=1:numel(d_m_unshuf)
    subplot(1,2,i);
    
    d = d_m_unshuf(i).data;
    Z=linkage(d,'complete','cityblock');
    f=figure;
    [~, ~, p_rows]=dendrogram(Z,0);
    close(f);
    Z=linkage(d','complete','cityblock');
    f=figure;
    [~, ~, p_cols]=dendrogram(Z,0);
    close(f);
    
    
    imagesc(d(p_rows,p_cols)');
    caxis([0 1]);
    colormap(flip(bone));
    cb = colorbar;
    cb.Label.String = 'Bootstrap probability';
    title(sprintf('decathlon batch #%i',i));
    set(gca,'XTick',1:numel(d_m_unshuf(i).rows),'XTickLabels',...
        d_m_unshuf(i).rows(p_rows),'XTickLabelRotation',90,...
        'YTick',1:numel(d_m_unshuf(i).cols),'YTickLabels',...
        pretty_labels(d_m_unshuf(i).cols(p_cols)),'FontSize',8,'TickLength',[0 0]);
end


% plot matrices
figure;
for i=1:numel(d_m_shuf)
    subplot(1,2,i);
    
    d = d_m_shuf(i).data;
    Z=linkage(d,'average','correlation');
    f=figure;
    [~, ~, p_rows]=dendrogram(Z,0);
    close(f);
    Z=linkage(d','average','correlation');
    f=figure;
    [~, ~, p_cols]=dendrogram(Z,0);
    close(f);
    
    imagesc(d(p_rows,p_cols)');
    colorbar;
    set(gca,'XTick',1:numel(d_m_unshuf(i).rows),'XTickLabels',...
        d_m_shuf(i).rows(p_rows),'XTickLabelRotation',90,...
        'YTick',1:numel(d_m_unshuf(i).cols),'YTickLabels',...
        pretty_labels(d_m_unshuf(i).cols(p_cols)),'FontSize',8,'TickLength',[0 0]);
end


%% read in data

fpaths = recursiveSearch(fdir,'ext','.csv','keyword','d');

[dd,ff,ee] = cellfun(@(s) fileparts(s), fpaths, 'UniformOutput', false);
batch_idx = regexp(ff,'(?<=d)[0-9]','match');
batch_idx = cat(1,batch_idx{:});
batch_idx = cellfun(@str2double,batch_idx);
grp_names = regexp(ff,'(?<=unshuffled_).*(?=_cat_x_gene)','match');
filt = cellfun(@isempty,grp_names);
batch_idx(filt) = [];
fpaths(filt) = [];
grp_names(filt) = [];

all_metric_idx = find(strcmp(cat(1,grp_names{:}),'all'));
batch_idx(all_metric_idx) = [];
bg_path = fpaths(all_metric_idx);
grp_names(all_metric_idx) = [];
batch_names = arrayfun(@(i) grp_names(batch_idx==i), 1:numel(bg_path), 'UniformOutput', false);
fpaths(all_metric_idx) = [];
batch_paths = arrayfun(@(i) fpaths(batch_idx==i), 1:numel(bg_path), 'UniformOutput', false);


%% sort data to common D1-D2 reference frame

% read in background categories
cats = cell(numel(bg_path),1);
genes = cell(numel(bg_path),1);
data = cell(numel(bg_path),1);
for j = 1:numel(bg_path)
    d_bg = read_csv_mat(bg_path{j});
    data{j} = d_bg.data./max(d_bg.data(:));
    data{j} = data{j}(:,1:end-1);
    cats{j} = d_bg.rows;
    genes{j} = d_bg.cols(1:end-1);
end

for i=1:numel(cats)
    cats{i} = cellfun(@(s) lower(s), cats{i}, 'UniformOutput', false);
end
for k=1:numel(cats)
    for i=1:numel(cats{k})
        tmp_cat = cats{k}{i};
        for j=1:numel(key_find)
            start_idx = strfind(tmp_cat,key_find{j});
            if start_idx
                tmp_cat(start_idx:start_idx+numel(key_find{j})-1) = '';
                tmp_cat = cat(2,tmp_cat(1:start_idx-1),key_replace{j},tmp_cat(start_idx:end));
            end
        end
        cats{k}{i} = strtrim(tmp_cat);
    end
end



figure;
gene_mask = cell(numel(data),1);
row_perm = cell(numel(data),1);
col_perm = cell(numel(data),1);
for i = 1:numel(data)
    
    data{i}(isnan(data{i})) = 0;
    [~,ii] = max(data{i});
    gene_mask{i} = false(size(data{i}));
    gene_mask{i}(sub2ind(size(data{i}),ii,1:size(data{i},2))) = true;
    
    % sort rows and cols
%   row_avg = cellfun(@(dd) nanmean(dd), num2cell(data{i},2), 'UniformOutput', false);
%   row_avg = cat(1,row_avg{:});
%   [~,p1] = sort(row_avg,'descend');
%   [~,p1] = sort(sum(data{i}>0.05,2),'descend')
    [~,row_perm{i}] = ismember(flip(d_p_unshuf(i).rows),cats{i});
    [r,c] = find(gene_mask{i}(row_perm{i},:));
    [~,col_perm{i}] = sort(r);
    col_perm{i} = [c(col_perm{i}); c(~ismember(c,c(col_perm{i})))];
    cats{i} = cats{i}(row_perm{i});
    genes{i} = genes{i}(col_perm{i});
    data{i} = data{i}(row_perm{i},col_perm{i});
    
    subplot(numel(data),1,i);
    imagesc(data{i});
    caxis([0 1]);
    colormap(flip(bone,1));
    set(gca,'YTick',1:numel(cats{i}),'YTickLabel',cats{i},'TickLength',[0 0]);
    ax = gca;
    ah.XAxis.FontSize = 6;
    ax.YAxis.FontSize = 5;
    ax.YAxis.Label.FontSize = 8;
    cb = colorbar;
    cb.Label.FontSize = 6;
    cb.FontSize = 6;
end



%% plot heatmaps partitioned by a priori groups
    

for j = 1:numel(bg_path)
    f = figure;
    for i = 1:numel(batch_paths{1})
        subplot(ceil(numel(batch_paths{1})/2),2,i);
        d_bg = read_csv_mat(bg_path{j});
        d_p = read_csv_mat(batch_paths{j}{i});
        d_p.data = d_p.data./max(d_p.data(:));
        imagesc(d_p.data);
        colormap(flip(bone));
        colorbar;
        caxis([0 1]);
        title(sprintf('Decathlon %i - %s',j,batch_names{j}{i}{1}));
        ylabel('KEGG enrichment');
        xlabel('gene');
        set(gca,'YTick',1:numel(cats{j}),'YTickLabel',cats{j},'TickLength',[0 0],'XTick',[]);
        ax = gca;
        ax.YAxis.FontSize = 5;
        ax.YAxis.Label.FontSize = 10;
    end   
end

%%

gene_hits.fbgn = cell(numel(has_match),1);
gene_hits.kegg = cell(numel(has_match),1);
gene_hits.name = cell(numel(has_match),1);

gene_hits.fbgn(has_match) = fbgn;
gene_hits.kegg = a;
gene_hits.name(has_match) = fb_gene_name;

%%

idx = 1;
d = num2cell(data{idx} > 0,1);
gene_lists = cellfun(@(ii) gene_hits(idx).name(ii), d, 'UniformOutput', false);

cat_genes = cell(numel(cats),1);
gnames = gene_hits(idx).name;
for i=1:numel(cats{idx})
    d = data{idx}(i,:);
    tmp = gnames(d>0);
    tmp = cellfun(@(s) sprintf('%s, ',s),tmp,'UniformOutput',false);
    [~,ii] = sort(lower(tmp));
    tmp = tmp(ii);
    tmp = cat(2,tmp{:});
    tmp(end-1:end) = [];
    cat_genes{i} = tmp;
    fprintf('%s\n',cats{idx}{i});
    fprintf('%s\n\n',cat_genes{i});
end

%%

D_enrichment_results = repmat(struct(...
    'experiment',[],'isogenic',[],'prob_gene_given_cat',[],...
    'cat_labels',[],'gene_labels',[],'gene_kegg_ids',[],'gene_fb_gene_num',[],...
    'avg_min_pval',[],'avg_metrics_hit',[]),numel(data),1);
exps = arrayfun(@(i) sprintf('decathlon-%i',i), 1:numel(data), 'UniformOutput', false);
is_iso = true(numel(data),1);

for i=1:numel(D_enrichment_results)
    D_enrichment_results(i).experiment = exps{i};
    D_enrichment_results(i).isogenic = is_iso(i);
    D_enrichment_results(i).prob_gene_given_cat = data{i};
    D_enrichment_results(i).cat_labels = cats{i};
    D_enrichment_results(i).gene_labels = gene_hits(i).name;
    D_enrichment_results(i).gene_kegg_ids = gene_hits(i).kegg;
    D_enrichment_results(i).gene_fb_gene_num = gene_hits(i).fbgn;
    D_enrichment_results(i).avg_min_pval = flip(d_p_unshuf(i).data);
    D_enrichment_results(i).avg_metrics_hit = flip(d_n_unshuf(i).data);
end

