
D_p = pair_rnaseq(D_seq(1:2));
pairs = unique_idx_pairs(numel(D_p),1);
figure;
for i=1:size(pairs,1)
    
    % get expression data and filter out empty reads
    d1 = D_p(pairs(i,1)).data;
    filt = nansum(d1,2)<0.5E6;
    d1 = d1(~filt,:);
    d2 = D_p(pairs(i,2)).data;
    filt = nansum(d2,2)<0.5E6;
    d2 = d2(~filt,:);
    
    % normalize to reads per million
    tot_d1 = sum(d1,2);
    scale_d1 = tot_d1./min(tot_d1);
    d1 = d1./repmat(scale_d1,1,size(d1,2));
    tot_d2 = sum(d2,2);
    scale_d2 = tot_d2./min(tot_d2);
    d2 = d2./repmat(scale_d2,1,size(d2,2));
    
    % compute mean and variance of reads for each gene
    gene_var = cellfun(@var,{d1;d2},'UniformOutput',false);
    gene_mu = cellfun(@mean,{d1;d2},'UniformOutput',false);

    % plot distribution of variance in each dataset
    subplot(2,3,(i-1)*3+1); hold on;
    [kde,x] = cellfun(@(v) ksdensity(log10(v),linspace(-10,10,50)),...
        gene_var, 'UniformOutput', false);
    lh = gobjects(2,1);
    lh(1) = plot(x{1},kde{1},'Color',[.1 .1 0.8],'LineWidth',2);
    lh(2) = plot(x{2},kde{2},'Color',[1 .5 0],'LineWidth',2);
    xlabel('gene expression variance');
    ylabel('PDF');
    leg_labels = arrayfun(@(i) sprintf('batch %i',i),...
        pairs(i,:), 'UniformOutput', false);
    legend(lh,leg_labels,'Location','NorthWest');
    title('Raw');
    
    % variance scatter plot
    subplot(2,3,(i-1)*3+2); hold on;
    pretty_scatter(log10(gene_mu{1}),log10(gene_mu{2}),'k','MarkerSize',1);
    xlabel(sprintf('variance (batch %i)',pairs(i,1)));
    ylabel(sprintf('variance (batch %i)',pairs(i,2)));
    title('Raw');

    % mean-variance plot
    subplot(2,3,(i-1)*3+3); hold on;
    lh = gobjects(2,1);
    lh(1) = pretty_scatter(log10(gene_mu{1}),log10(gene_var{1}),[.1 .1 0.8],'MarkerSize',1);
    lh(2) = pretty_scatter(log10(gene_mu{2}),log10(gene_var{2}),[1 .5 0],'MarkerSize',1);
    xlabel({'gene expression';'log(mean)'});
    ylabel({'gene expression';'log(variance)'});
    set(gca,'XLim',[-2 6],'YLim',[-2 12]);
    title('Raw');
    
    % get expression data and filter out empty reads
    d1 = D_p(pairs(i,1)).data;
    filt = nansum(d1,2)<0.5E6;
    d1 = d1(~filt,:);
    d2 = D_p(pairs(i,2)).data;
    filt = nansum(d2,2)<0.5E6;
    d2 = d2(~filt,:);
    
    % normalize to reads per million
    tot_d1 = sum(d1,2);
    scale_d1 = tot_d1./min(tot_d1);
    d1 = d1./repmat(scale_d1,1,size(d1,2));
    tot_d2 = sum(d2,2);
    scale_d2 = tot_d2./min(tot_d2);
    d2 = d2./repmat(scale_d2,1,size(d2,2));
    
    % quantile normalize
%     d1 = quantile_normalize(d1')';
%     d2 = quantile_normalize(d2')';
    q_norm = quantile_normalize([d1;d2]')';
    gene_filt = nanmean(q_norm)>0;
    d1 = q_norm(1:size(d1,1),gene_filt);
    d2 = q_norm(size(d1,1)+1:end,gene_filt);    
    
    % compute mean and variance of reads for each gene
    gene_var = cellfun(@var,{d1;d2},'UniformOutput',false);
    gene_mu = cellfun(@mean,{d1;d2},'UniformOutput',false);

    % plot distribution of variance in each dataset
    subplot(2,3,i*3+1); hold on;
    [kde,x] = cellfun(@(v) ksdensity(log10(v),linspace(-10,10,50)),...
        gene_var, 'UniformOutput', false);
    lh = gobjects(2,1);
    lh(1) = plot(x{1},kde{1},'Color',[.1 .1 0.8],'LineWidth',2);
    lh(2) = plot(x{2},kde{2},'Color',[1 .5 0],'LineWidth',2);
    xlabel('gene expression variance');
    ylabel('PDF');
    leg_labels = arrayfun(@(i) sprintf('batch %i',i),...
        pairs(i,:), 'UniformOutput', false);
    legend(lh,leg_labels,'Location','NorthWest');
    title('Quantile Normalized');
    
    % variance scatter plot
    subplot(2,3,i*3+2); hold on;
    pretty_scatter(log10(gene_var{1}),log10(gene_var{2}),'k','MarkerSize',1);
    xlabel(sprintf('variance (batch %i)',pairs(i,1)));
    ylabel(sprintf('variance (batch %i)',pairs(i,2)));
    title('Quantile Normalized');


    % mean-variance plot
    subplot(2,3,i*3+3); hold on;
    lh = gobjects(2,1);
    lh(1) = pretty_scatter(log10(gene_mu{1}),log10(gene_var{1}),[.1 .1 0.8],'MarkerSize',1);
    lh(2) = pretty_scatter(log10(gene_mu{2}),log10(gene_var{2}),[1 .5 0],'MarkerSize',1);
    xlabel({'gene expression';'log(mean)'});
    ylabel({'gene expression';'log(variance)'});
    set(gca,'XLim',[-2 6],'YLim',[-2 12]);
    title('Quantile Normalized');
end

%%
figure;
for i=1:numel(D_p)
    subplot(2,numel(D_p),i);
    
    % get expression data and filter out empty reads
    d = D_p(i).data;
    filt = nansum(d,2)<0.2E6;
    d = d(~filt,:);
    cts = nanmean(d);
    [kde,x] = ksdensity(log10(cts),linspace(-5,5,100));
    plot(x,kde,'k-','Linewidth',2);
    xlabel('log(mean reads per gene)');
    ylabel('PDF');
    title(sprintf('Unnormalized (batch %i)',i));
end

% get expression data and filter out empty reads
d1 = D_p(1).data;
filt = nansum(d1,2)<0.5E6;
d1 = d1(~filt,:);
d2 = D_p(2).data;
filt = nansum(d2,2)<0.5E6;
d2 = d2(~filt,:);
q_norm = quantile_normalize([d1;d2]')';
d1 = q_norm(1:size(d1,1),:);
d2 = q_norm(size(d1,1):end,:);



subplot(2,numel(D_p),3);
cts = nanmean(d1);
[kde,x] = ksdensity(log10(cts),linspace(-5,5,100));
plot(x,kde,'k-','Linewidth',2);
xlabel('log(mean reads per gene)');
ylabel('PDF');
title(sprintf('Quantile Normalized (batch %i)',1));

subplot(2,numel(D_p),4);
cts = nanmean(d2);
[kde,x] = ksdensity(log10(cts),linspace(-5,5,100));
plot(x,kde,'k-','Linewidth',2);
xlabel('log(mean reads per gene)');
ylabel('PDF');
title(sprintf('Quantile Normalized (batch %i)',2));

%%

% load behavioral data
D_b = load_decathlon_structs(fdir,'D123_als');
D_b = D_b(1:2);
D_b = pair_decathlon_structs(D_b,'CollapseMode','PCA','CollapseFields','none');
f = true(192,1);
f([99 114])=false;
d = NaN(size(D_b(1).data));
d(f,:) = D_b(1).data;
D_b(1).data = d;

% pair, filter, and quantile normalize expression data
D_p = pair_rnaseq(D_seq(1:2));
for i=1:numel(D_p)
    d = D_p(1).data;
    fly_filt = nansum(d,2)>0.5E6;
    d = d(fly_filt,:);
    
    % quantile normalize
    q_norm = quantile_normalize(d')';
    gene_filt = nanmean(q_norm)>10;
    
    % update data
    D_p(i).data = q_norm(:,gene_filt);
    D_p(i).ID = D_p(i).ID(fly_filt) - (i-1)*192;
    D_p(i).geneID = D_p(i).geneID(gene_filt);
    D_b(i).data = D_b(i).data(D_p(i).ID,:);
end

% train gene x behavior linear models
model_r_squared = cell(numel(D_p),1);
model_p_value = cell(numel(D_p),1);
shuf_r_squared = cell(numel(D_p),1);
shuf_p_value = cell(numel(D_p),1);
for i=1:numel(D_p)
    
    model_r_squared{i} = NaN(size(D_p(i).data,2),size(D_b(i).data,2));
    model_p_value{i} = NaN(size(D_p(i).data,2),size(D_b(i).data,2));
    shuf_r_squared{i} = NaN(size(D_p(i).data,2),size(D_b(i).data,2));
    shuf_p_value{i} = NaN(size(D_p(i).data,2),size(D_b(i).data,2));
    shuf_b = shuffle_columns(D_b(i).data);
    
    
    for j=1:size(D_b(i).data,2)
        % print notification
        fprintf('D%i, Behavior %i of %i\n',i,j,size(D_b(i).data,2));

        obs_r = NaN(10,size(D_p(i).data,2));
        shuf_r = NaN(10,size(D_p(i).data,2));
        for k = 1:size(obs_r,1)
            rnd_idx = randperm(size(D_b(i).data,1),ceil(0.8*size(D_b(i).data,1)));
            obs_r(k,:) = corr(D_b(i).data(rnd_idx,j),D_p(i).data(rnd_idx,:));
            shuf_r(k,:) = corr(shuf_b(rnd_idx,j),D_p(i).data(rnd_idx,:));
        end
        obs_r = mean(obs_r);
        shuf_r = mean(shuf_r);
        
        [~,model_p_value{i}(:,j)] = corr(D_b(i).data(:,j),D_p(i).data);
        [~,shuf_p_value{i}(:,j)] = corr(shuf_b(:,j),D_p(i).data); 
    end
end

%% Create a FBgno. to KEGG ID lookup table

% export background
chunk_sz = 4000;
all_genes = D_seq(1).geneID;
for i=1:ceil(numel(all_genes)/chunk_sz)
    
    if i==ceil(numel(all_genes)/chunk_sz)
        write_genes = all_genes((i-1)*chunk_sz+1:end);
    else
        write_genes = all_genes((i-1)*chunk_sz+1:chunk_sz*i);
    end
    fname = sprintf('%sbackground_%i.csv',save_dir,i);
    fid = fopen(fname,'W+');
    cellfun(@(s) fprintf(fid,'%s\n',s), write_genes);
    fclose(fid);
end


%% Export significant hits for each apriori group

% load kegg gene id lookup table
lookup_path = ['C:\Users\winsl0w\Documents\decathlon\decathlon_paper_data\'...
    'decathlon_go_geneid_data\apriori_grp_hits\background\'...
    'geneid_lookup_table.mat'];
load(lookup_path);

% set significance threshold
alpha = 0.05;
save_dir = ['C:\Users\winsl0w\Documents\decathlon\decathlon_paper_data\'...
    'decathlon_go_geneid_data\apriori_grp_hits\kegg_gid\'];

% export background to file
fname = sprintf('%sbackground_kegg_gid.csv',save_dir);
fid = fopen(fname,'W+');
bg_kegg_ids = fbgn2kegg(D_p(1).geneID,geneid_lookup_tbl);
cellfun(@(s) fprintf(fid,'%s\n',s), bg_kegg_ids);
fclose(fid);

% query field names and group aprior fields
f = D_b(1).fields;
[~, grp_name, grp_idx] = group_apriori_fields(D_b(1));

% iterate over decathlon batches and aprior groups
for i=1:numel(D_b)
    for j=1:numel(grp_name)
        % select gene hits
        grp_hits = model_p_value{i}(:,grp_idx{j}) < alpha;
        grp_genes = D_p(i).geneID(any(grp_hits,2));
        
        % convert flybase gene no. to kegg gene ids
        kegg_ids = fbgn2kegg(grp_genes,geneid_lookup_tbl);
        
        % export hits to file
        grp_name{j}(grp_name{j}==' ') = '_';
        fname = sprintf('%sD%i_%s_kegg_gid.csv',save_dir,i,grp_name{j});
        fid = fopen(fname,'W+');
        cellfun(@(s) fprintf(fid,'%s\n',s), kegg_ids);
        fclose(fid);
    end
end

%% Export significant hits for each apriori group

% load kegg gene id lookup table
lookup_path = ['C:\Users\winsl0w\Documents\decathlon\decathlon_paper_data\'...
    'decathlon_go_geneid_data\apriori_grp_hits\background\'...
    'geneid_lookup_table.mat'];
load(lookup_path);

% set significance threshold
alpha = 0.05;
save_dir = ['C:\Users\winsl0w\Documents\decathlon\decathlon_paper_data\'...
    'decathlon_go_geneid_data\apriori_grp_hits\kegg_gid\'];

% export background to file
fname = sprintf('%sbackground_kegg_gid.csv',save_dir);
fid = fopen(fname,'W+');
bg_kegg_ids = fbgn2kegg(D_p(1).geneID,geneid_lookup_tbl);
cellfun(@(s) fprintf(fid,'%s\n',s), bg_kegg_ids);
fclose(fid);

% query field names and group aprior fields
f = D_b(1).fields;
[~, grp_name, grp_idx] = group_apriori_fields(D_b(1));

% iterate over decathlon batches and aprior groups
for i=1:numel(D_b)
    for j=1:numel(f)
        % select gene hits
        grp_hits = model_p_value{i}(:,j) < alpha;
        grp_genes = D_p(i).geneID(any(grp_hits,2));
        
        % convert flybase gene no. to kegg gene ids
        kegg_ids = fbgn2kegg(grp_genes,geneid_lookup_tbl);
        
        % export hits to file
        f{j}(f{j}==' ') = '_';
        fname = sprintf('%sD%i_%s_kegg_gid.csv',save_dir,i,f{j});
        fid = fopen(fname,'W+');
        cellfun(@(s) fprintf(fid,'%s\n',s), kegg_ids);
        fclose(fid);
    end
end


%% plot distributions

figure;
n_plots = 5;
for i=1:numel(D_p)
    subplot(n_plots,numel(D_p),i); hold on;
    bins = linspace(0,1,100);
    ct = histc(model_p_value{i}(:),bins);
    lh_obs = plot(bins,ct,'r-','Linewidth',1.5);
    ct = histc(shuf_p_value{i}(:),bins);
    lh_shuf = plot(bins,ct,'k--','Linewidth',1.5);
    xlabel('-log(p-value)');
    ylabel('PDF');
    title(sprintf('Decathlon-%i',i));
    
    subplot(n_plots,numel(D_p),i+2); hold on;
    imagesc(-log10(model_p_value{i})');
    cb = colorbar;
    caxis([0 4]);
    axis tight
    xlabel('genes');
    ylabel('behaviors');
    cb.Label.String = '-log(p-value)';
    
    subplot(n_plots,numel(D_p),i+4); hold on;
    imagesc(-log10(shuf_p_value{i})');
    cb = colorbar;
    caxis([0 4]);
    axis tight
    xlabel('genes');
    ylabel('behaviors');
    cb.Label.String = '-log(p-value)';
    
    subplot(n_plots,numel(D_p),i+6); hold on;
    n_sig_obs = sum(model_p_value{i}<0.05)./size(model_p_value{i},1);
    n_sig_shuf = sum(shuf_p_value{i}<0.05)./size(model_p_value{i},1);
    [kde,x] = ksdensity(n_sig_obs,linspace(0,0.25,100));
    plot(x,kde,'r');
    [kde,x] = ksdensity(n_sig_shuf,linspace(0,0.25,100));
    plot(x,kde,'k');
    xlabel('fraction significant genes per behavior');
    ylabel('PDF');
    
    subplot(n_plots,numel(D_p),i+8); hold on;
    n_sig_obs = sum(model_p_value{i}<0.05,2)./size(model_p_value{i},2);
    n_sig_shuf = sum(shuf_p_value{i}<0.05,2)./size(shuf_p_value{i},2);
    [kde,x] = ksdensity(n_sig_obs,linspace(0,0.2,20));
    plot(x,kde,'r');
    [kde,x] = ksdensity(n_sig_shuf,linspace(0,0.2,20));
    plot(x,kde,'k');
    xlabel('fraction significant behaviors per gene');
    ylabel('PDF');
end


