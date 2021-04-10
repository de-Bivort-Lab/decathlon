% Plot panels from the boostrapped KEGG enrichment analysis. 
% 1) Create bootstrapped shuffled and unshuffled average minimum <i>p</i>-value bar plots. 
% 2) Create bootstrapped average number metrics hit by average minimum <i>p</i>-value scatter plots.
% 3) Create average minimum <i>p</i>-value paired dot plots. 
% 4) Create KEGG pathway by gene boostrapped heatmaps split by a priori groups.

D = D_enrichment_results_us;

% generate p-value bar plots for shuffled and unshuffled data
figure;
for i=1:numel(D)
    subplot(numel(D),1,i);  hold on;
    h1 = barh(flip(D(i).avg_min_pval),'FaceColor',[.6 .6 .6], 'EdgeColor', 'none');
    h2 = barh(flip(D(i).shuffled_avg_min_pval),'FaceColor',[0 0 0], 'EdgeColor', 'none');   
    set(gca,'YTick',1:numel(D(i).cat_labels),...
        'YTickLabel',flip(D(i).cat_labels),'TickLength',[0 0],'XLim',[0 15]);
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
for i=1:numel(D)
    
    % scatter unshuffled data
    subplot(numel(D),1,i);  hold on;
    lh1 = pretty_scatter(D(i).avg_metrics_hit,D(i).avg_min_pval,[.6 .6 .6],'MarkerSize',3);
    title(sprintf('decathlon %i',i));
    xlabel('num metrics');
    ylabel('-log_{10}[p]');
    
    % apply text labels to significant nodes
    p_filt = find(D(i).avg_min_pval>2);
    x = D(i).avg_metrics_hit;
    y = D(i).avg_min_pval;
    for j=1:numel(p_filt)
        text(x(p_filt(j))-.2,y(p_filt(j)),D(i).cat_labels(p_filt(j)),...
            'HorizontalAlignment','right','FontSize',6);
    end
end


% pair dot plot
pairs = unique_idx_pairs(numel(D),true);
D_p = D;
for i=1:size(pairs,1)
    [p_a,p_b] = ismember(D_p(pairs(i,1)).cat_labels,D_p(pairs(i,2)).cat_labels);
    D_p(pairs(i,1)).avg_min_pval = D_p(pairs(i,1)).avg_min_pval(p_a);
    D_p(pairs(i,2)).avg_min_pval = D_p(pairs(i,2)).avg_min_pval(p_b(p_b>0));
    D_p(pairs(i,1)).cat_labels = D_p(pairs(i,1)).cat_labels(p_a);
    D_p(pairs(i,2)).cat_labels = D_p(pairs(i,2)).cat_labels(p_b(p_b>0));
end
figure;
for i=1:size(pairs,1)
    subplot(1,size(pairs,1),i);
    a = D_p(pairs(i,1)).avg_min_pval;
    b = D_p(pairs(i,2)).avg_min_pval;
    plot_pair_dot(a,b,'k',[]);
    xticks = arrayfun(@(ii) sprintf('D%i',ii),pairs(i,:),'UniformOutput',false);
    set(gca,'XTick',[1 2],'XTickLabel',xticks,'YLim',[0 15]);
    title('unshuffled');
    ylabel('-log_{10}[p-value]');
    hold on;
    for j=1:numel(D_p(pairs(i,2)).cat_labels)
       text(2.1,b(j),D_p(pairs(i,2)).cat_labels{j},'FontSize',6); 
    end
end

%%

% ----  PLOT KEGG PATHWAY x BEHAVIORAL METRIC HEATMAPS ---- %
% load behavioral data
D_b = load_decathlon_structs(pwd,'D_us');

% get field sorting order
for i=1:numel(D_b)
    % sort by apriori group for d1, then match d2 to d1
    [~,grp_idx,fields] = group_unsupervised_clusters;
    D_b(i).fields = fields;
    D_b(i).data = D_b(i).pdfs(:,cat(2,grp_idx{:}));
    
    [~, ~, grp_idx] = group_apriori_fields(D_b(i));
    
    figure;
    imagesc(D(i).p_metric);
    colormap(flip(bone));
    caxis([0 1]);
    cb = colorbar;
    cb.Label.String = 'Bootstrap probability';
    title(sprintf('decathlon batch #%i',i));
    set(gca,'YTick',1:numel(D(i).cat_labels),'YTickLabels',...
        D(i).cat_labels,'FontSize',8,'TickLength',[0 0],...
        'XTick',1:numel(D(i).metric_labels),'XTickLabels',...
        pretty_labels(D(i).metric_labels),...
        'XTickLabelRotation',90);
    
    id_group = cellfun(@(gi,i) repmat(i,numel(gi),1),...
        grp_idx, num2cell(1:numel(grp_idx))', 'UniformOutput', false);
    id_group = cat(1,id_group{:});
    
    % create patch coordinates for groups
    % create patches for a priori groups and assays
    vx_groups = repmat([0;0;1;1;0],1,numel(grp_idx));
    vy_groups = NaN(5,numel(grp_idx));
    color_groups = jet(numel(grp_idx));
    for j = 1:numel(grp_idx)
        idx = [find(id_group==j,1) find(id_group==j,1,'Last')];
        if ~isempty(idx)
            y = [idx(1) idx([2 2])+1 idx([1 1])]';
            y = y + 0.2.*[1;-1;-1;1;1];
            vy_groups(:,j) = y;
        end
    end

    max_y = numel(cat(1,grp_idx{:}));
    hold on;
    patch('YData',vx_groups-1,'XData',vy_groups-0.5,...
        'FaceColor','flat','CData',linspace(0,1,numel(grp_idx)));
    axis(gca,'equal','tight');
end

% ---- PLOT KEGG PATHWAY BY GENE HEATMAPS ---- %
figure;
for i = 1:numel(D)

    subplot(numel(D),1,i);
    imagesc(D(i).prob_gene_given_cat);
    caxis([0 1]);
    colormap(flip(bone,1));
    set(gca,'YTick',1:numel(D(i).cat_labels),'YTickLabel',D(i).cat_labels,'TickLength',[0 0]);
    ax = gca;
    ah.XAxis.FontSize = 6;
    ax.YAxis.FontSize = 5;
    ax.YAxis.Label.FontSize = 8;
    cb = colorbar;
    cb.Label.FontSize = 6;
    cb.FontSize = 6;
end

% ---- PLOT KEGG PATHWAY BY GENE HEATMAPS PARTITIONED BY APRIORI GRP  ---- %
for j = 1:numel(D)
    f = figure;
    fn = fieldnames(D(j).apriori_cat_x_gene);
    for i = 1:numel(fn)
        subplot(ceil(numel(fn)/2),2,i);
        imagesc(D(j).apriori_cat_x_gene.(fn{i}));
        colormap(flip(bone));
        colorbar;
        caxis([0 1]);
        title(sprintf('Decathlon %i - %s',j,fn{i}));
        ylabel('KEGG enrichment');
        xlabel('gene');
        set(gca,'YTick',1:numel(D(j).cat_labels),'YTickLabel',D(j).cat_labels,'TickLength',[0 0],'XTick',[]);
        ax = gca;
        ax.YAxis.FontSize = 5;
        ax.YAxis.Label.FontSize = 10;
    end   
end

