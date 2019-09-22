
for j=1:4
    
    ah = subplot(2,2,j);
    switch j
        case 1
            collapse_mode = 'average';
            collapse_fields = 'none';
            title_str = 'circ expanded';
        case 2
            collapse_mode = 'average';
            collapse_fields = 'circadian';
            title_str = 'circ collapsed';
        case 3
            collapse_mode = 'average';
            collapse_fields = 'all';
            title_str = 'a priori collapsed';
        case 4
            collapse_mode = 'PCA';
            collapse_fields = 'all';
            title_str = 'a priori PCA';
    end
    D_col = D;
    for i=1:numel(D_col)
        D_col(i) = collapseMetrics(D_col(i),'CollapseMode',collapse_mode,...
            'CollapseFields',collapse_fields,'PCs',2);
    end
    colors = [0 0 0.5; 0.5 0 0; 0 0.5 0];
    for i=1:numel(D_col)
        pair_idx = upper_triangle_idx(numel(D_col(i).fields));
        [~,p_tmp] = corr(D_col(i).data,'Type','spearman','rows','pairwise');
        hold on;
        ph = histogram(ah,p_tmp(pair_idx),linspace(0,1,20),...
            'FaceColor',colors(i,:),'FaceAlpha',0.5);
        uistack(ph,'bottom');
        xlabel('p-value');
        set(gca,'XLim',[0 1]);
    end
    hold on;
    %plot([0 1],[0 1],'k--','LineWidth',0.5);

    title(sprintf('correlation matrix p-values (%s)',title_str));
    labels = arrayfun(@(i) sprintf('D%i',i), numel(D_col):-1:1, 'UniformOutput', false);
    legend([labels],'location','NorthEast');
end