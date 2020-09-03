% Plot heatmaps of decathlon RNAseq data

% set params
min_tot_reads = 1E5;
min_rpm = 10;
molaspass=interp1([1 51 102 153 204 256],...
        [0 0 0; 0 0 .75; .5 0 .8; 1 .1 0; 1 .9 0; 1 1 1],1:256);

% filter and sort genes to align decathlon-1 and decathlon-2 data
D_p = pair_rnaseq(D_seq(1:2));

% plot heatmaps
figure;
reads = cell(numel(D_p),1);
for i=1:numel(reads)
    reads{i} = D_p(i).data;
    
    % filter out flies with too few reads
    fly_filt = nansum(reads{i},2) > min_tot_reads;
    reads{i} = reads{i}(fly_filt,:);
    
    % normalize to reads per million
    tot_reads = sum(reads{i},2);
    scale_factor = tot_reads ./ min(tot_reads);
    reads{i} = reads{i} ./ repmat(scale_factor,1,size(reads{i},2));
    
    % plot bootstrap pca
    subplot(1,numel(reads),i);
    plot_pca_bootstrap(reads{i},50);
    set(gca,'XLim',[1 100],'XScale','log');
end

% quantile normalize reads
all_reads = cat(1,reads{:});
gene_filt = nanmean(all_reads) >= min_rpm;
filt_idx = sum(~gene_filt);
q_norm = quantile_normalize(cat(1,reads{:})')';

% parse reads and plot heatmaps
figure;
q_reads = cell(numel(reads),1);
for i=1:numel(q_reads)
    n_rows = [0;cellfun(@(x) size(x,1), reads(1:i))];
    q_reads{i} = q_norm(sum(n_rows(1:i))+1:sum(n_rows),:);


    % find row and column permutation to sort by expression level
    [~,perm_cols] = sort(nanmean(q_reads{i}));
    [~,perm_rows] = sort(nanmean(q_reads{i},2));
    x = sum(nanmean(q_reads{i})<10);
    
    
    % plot
    subplot(2,1,i); hold on;
    imagesc(log10(q_reads{i}(perm_rows,perm_cols)));
    plot([x x],[1 size(q_reads{i},1)],'--','Color',[.7 .7 .7],'LineWidth',1);
    colormap(molaspass);
    cb = colorbar;
    caxis([0 2]);
    cb.Label.String = 'log[RPM]';
    xlabel('genes');
    ylabel('individual flies');
    title(sprintf('decathlon %i',i));
    set(gca,'XLim',[1 size(reads{i},2)]);
end


%%

figure;
for i=1:numel(q_reads)
    % plot bootstrap pca
    subplot(1,numel(reads),i);
    plot_pca_bootstrap(q_reads{i},50);
    set(gca,'XLim',[1 100],'XScale','log');
end