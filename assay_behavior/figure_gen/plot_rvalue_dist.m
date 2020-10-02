function [rvals,avg_pval_kde,ah,pvals] = plot_rvalue_dist(D,ah)

rvals = cell(size(D,2),1);
pvals = cell(size(D,2),1);
avg_pval_kde = cell(size(D,2),2);
r_ci95s = cell(size(D,2),1);
p_ci95s = cell(size(D,2),1);
hold(ah(1),'on');
hold(ah(2),'on');
colors = [.7 0 .7; 0 .7 0; .35 .35 .35];
nbins = 500;
lhs = gobjects(2,3);

n = cellfun(@(d) size(d,1), D(1,:));
for i=1:size(D,2)
    bs_data = D(:,i);
    [r,p] = cellfun(@(d) corr(d,'Type','spearman'),bs_data, 'UniformOutput', false);
    idx = cellfun(@(rr) upper_triangle_idx(size(rr,1)), r, 'UniformOutput', false);
    r = cellfun(@(rr,ii) rr(ii), r, idx, 'UniformOutput', false);
    p = cellfun(@(pp,ii) pp(ii), p, idx, 'UniformOutput', false);
    for j=1:numel(r)
        r{j}(r{j}==1) = 1-eps;
        p{j}(p{j}==1) = 1-eps;
    end

    bins = linspace(-1,1,nbins);
    kde = cellfun(@(rr) ksdensity(rr,bins,'Support',[-1 1]), r,...
        'UniformOutput', false);
    kde = cat(1,kde{:});
    avg_kde = median(kde);
    
    rvals{i} = r;
    pvals{i} = p;
    
    % compute and plot 95% CI
    ci95 = [prctile(kde,2.5);prctile(kde,97.5)];
    r_ci95s{i} = ci95;
    vx = [bins(1) bins bins(end) fliplr(bins)];
    vy = [ci95(2,1) ci95(1,:) ci95(2,end) fliplr(ci95(2,:))];
    patch('XData',vx(:),'YData',vy(:),'FaceColor',colors(i,:),'FaceAlpha',0.33,...
        'EdgeColor','none','Parent',ah(1));
    
    % plot distributions
    lhs(1,i) = plot(bins,avg_kde,'Color',colors(i,:),'Parent',ah(1),'LineWidth',1.5);

    % plot r-value KDE
    bins = linspace(0,1,nbins);
    p_kde = cellfun(@(pp) ksdensity(pp,bins,'Support',[-1 1]), p,...
        'UniformOutput', false);
    p_kde = cat(1,p_kde{:});
    avg_p_kde = mean(p_kde);
    avg_pval_kde{i,1} = avg_p_kde;
    
    % compute 95% CI
    ci95 = [prctile(p_kde,2.5);prctile(p_kde,97.5)];
    p_ci95s{i} = ci95;
    vx = [bins(1) bins bins(end) fliplr(bins)];
    vy = [ci95(2,1) ci95(1,:) ci95(2,end) fliplr(ci95(2,:))];
    patch('XData',vx(:),'YData',vy(:),'FaceColor',colors(i,:),'FaceAlpha',0.33,...
        'EdgeColor','none','Parent',ah(2));
    
    lhs(2,i) = plot(bins,avg_p_kde,'Color',colors(i,:),'Parent',ah(2),'LineWidth',1.5);
end

% compute for shuffled data sets
for i=size(D,2)+1:size(D,2)*2
    % calculate r-value and p-values
    bs_data = D(:,i-size(D,2));
    data = cellfun(@(d) shuffle_columns(d), bs_data, 'UniformOutput', false);
    [r,p] = cellfun(@(d) corr(d,'Type','spearman'),data, 'UniformOutput', false);
    idx = cellfun(@(rr) upper_triangle_idx(size(rr,1)), r, 'UniformOutput', false);
    r = cellfun(@(rr,ii) rr(ii), r, idx, 'UniformOutput', false);
    p = cellfun(@(pp,ii) pp(ii), p, idx, 'UniformOutput', false);
    for j=1:numel(r)
        r{j}(r{j}==1) = 1-eps;
        p{j}(p{j}==1) = 1-eps;
    end

    bins = linspace(-1,1,nbins);
    kde = cellfun(@(rr) ksdensity(rr,bins,'Support',[-1 1]), r,...
        'UniformOutput', false);
    kde = cat(1,kde{:});
    avg_kde = median(kde);

    rvals{i} = r;
    pvals{i} = p;
    
    % plot r-value KDE
    p_bins = linspace(0,1,nbins);
    p_kde = cellfun(@(pp) ksdensity(pp,p_bins,'Support',[-1 1]), p,'UniformOutput', false);
    p_kde = cat(1,p_kde{:});
    avg_p_kde = mean(p_kde);
    avg_pval_kde{i-size(D,2),2} = avg_p_kde;
    
    if i==3
    % compute and plot 95% CI
    ci95 = [prctile(kde,2.5);prctile(kde,97.5)];
    r_ci95s{i} = ci95;
    vx = [bins(1) bins bins(end) fliplr(bins)];
    vy = [ci95(2,1) ci95(1,:) ci95(2,end) fliplr(ci95(2,:))];
    patch('XData',vx(:),'YData',vy(:),'FaceColor',colors(i,:),'FaceAlpha',0.33,...
        'EdgeColor','none','Parent',ah(1));

    % plot distributions
    lhs(1,i) = plot(bins,avg_kde,'Color',colors(i,:),'Parent',ah(1),...
        'LineWidth',1.5,'LineStyle','--');
    xlabel(ah(1),'r-value');
    ylabel(ah(1),'PDF');
    legend(lhs(1,:),{'D12';'D3';'shuffled'});

    % compute 95% CI
    ci95 = [prctile(p_kde,2.5);prctile(p_kde,97.5)];
    vx = [p_bins(1) p_bins p_bins(end) fliplr(p_bins)];
    vy = [ci95(2,1) ci95(1,:) ci95(2,end) fliplr(ci95(2,:))];
    patch('XData',vx(:),'YData',vy(:),'FaceColor',colors(i,:),'FaceAlpha',0.33,...
        'EdgeColor','none','Parent',ah(2));

    lhs(2,i) = plot(p_bins,avg_p_kde,'Color',colors(i,:),'Parent',ah(2),...
        'LineWidth',1.5,'LineStyle','--');
    legend(lhs(2,:),{'D12';'D3';'shuffled'});
    xlabel(ah(2),'p-value');
    ylabel(ah(2),'PDF');
    end
end

% set axis limits
max_r = max(cellfun(@(c) max(c(:)), r_ci95s));
max_p = max(cellfun(@(c) max(c(:)), p_ci95s));
set(ah(1),'YLim',[0 max_r*1.1]);
set(ah(2),'YLim',[0 max_p*1.1]);







