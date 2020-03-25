
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
    
    % quantile normalize
    q_norm = quantile_normalize([d1;d2]')';
    gene_filt = nanmean(q_norm)>10;
    d1 = q_norm(1:size(d1,1),gene_filt);
    d2 = q_norm(size(d1,1):end,gene_filt);
%     d1 = quantile_normalize(d1')';
%     d2 = quantile_normalize(d2')';
    
    % compute mean and variance of reads for each gene
    gene_var = cellfun(@var,{d1;d2},'UniformOutput',false);
    gene_mu = cellfun(@mean,{d1;d2},'UniformOutput',false);

    % plot distribution of variance in each dataset
    subplot(size(pairs,1),3,(i-1)*3+1); hold on;
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
    
    % variance scatter plot
    subplot(size(pairs,1),3,(i-1)*3+2); hold on;
    pretty_scatter(log10(gene_var{1}),log10(gene_var{2}),'k');
    xlabel(sprintf('variance (batch %i)',pairs(i,1)));
    ylabel(sprintf('variance (batch %i)',pairs(i,2)));


    % mean-variance plot
    subplot(size(pairs,1),3,(i-1)*3+3); hold on;
    lh = gobjects(2,1);
    lh(1) = pretty_scatter(log10(gene_mu{1}),log10(gene_var{1}),[.1 .1 0.8],'MarkerSize',1);
    lh(2) = pretty_scatter(log10(gene_mu{2}),log10(gene_var{2}),[1 .5 0],'MarkerSize',1);
    xlabel({'gene expression';'log(mean)'});
    ylabel({'gene expression';'log(variance)'});
    set(gca,'XLim',[-2 6],'YLim',[-2 12]);
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


