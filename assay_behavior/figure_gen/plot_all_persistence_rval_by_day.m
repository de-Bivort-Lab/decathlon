%% group measures by assay
d_fields = standardize_fieldnames(d_fields);
[assay,metric,day] = parse_fieldnames(d_fields);
unique_assays = unique(assay);
unique_assays(strcmpi(unique_assays,'culling')) = [];
assay_grouped_data = cell(numel(unique_assays),1);
assay_grouped_fields = cell(numel(unique_assays),1);
for i=1:numel(unique_assays)
   assay_mask = strcmp(assay,unique_assays{i});
   assay_grouped_data{i} = data(assay_mask);
   assay_grouped_fields{i} = d_fields(assay_mask);
end

% iterate over each measure within each assay, group data by metric and
% compute the correlation matrix
data_by_metric = cell(numel(unique_assays),1);
labels_by_metric = cell(numel(unique_assays),1);
rval_by_metric = cell(numel(unique_assays),1);
rval_labels = cell(numel(unique_assays),1);
for i=1:numel(unique_assays)
   [~,metric,day] = parse_fieldnames(assay_grouped_fields{i});
   unique_metrics = unique(metric);
   rval_by_metric{i} = cell(numel(unique_metrics),1);
   data_by_metric{i} = cell(numel(unique_metrics),1);
   rval_labels{i} = unique_metrics;
   %metric_idx = str_list_contains(metric,unique_metrics);
   metric_idx = cellfun(@(s) strcmp(metric,s), unique_metrics, 'UniformOutput', false);
   metric_idx = cat(2,metric_idx{:});
   
   for j=1:numel(unique_metrics)
       metric_data = assay_grouped_data{i}(metric_idx(:,j));
       metric_labels = assay_grouped_fields{i}(metric_idx(:,j));
       metric_day = day(metric_idx(:,j));
       [~,day_sort] = sort(metric_day);
       metric_data = metric_data(day_sort);
       metric_labels = metric_labels(day_sort);
       data_by_metric{i}{j} = cat(2,metric_data{:});
       labels_by_metric{i}{j} = metric_labels;
       rval_by_metric{i}{j} = corr(data_by_metric{i}{j},'Type','spearman','rows','pairwise');
   end
end

%% Plot correlation matrices for each persistence assay

metric_permutation = {[6,5,2,1,4,3];[10,7,8,2,1,5,6,9,3,4];[7,4,5,2,1,6,3];...
    [8,5,6,2,1,7,3,4];[8,5,6,2,1,3,4,7]};
assay_labels = {'Circadian';'LED Y-Maze';'Optomotor';'Spatial Shade-light';'Y-Maze'};

for i =1:numel(data_by_metric)
    p = metric_permutation{i};
    d=data_by_metric{i}(p);
    L = labels_by_metric{i}(p);
    d = cat(2,d{:});
    L = cat(1,L{:});
    
    
    plotCorr(d,'Cluster',false,'Options',{'rows','pairwise'});
    step_sz = (size(d,2)/numel(rval_labels{i}));
    yticks = step_sz/2:step_sz:size(d,2);
    ah = gca;
    set(ah,'YTick',yticks,'YTickLabels',pretty_labels(rval_labels{i}(p)),...
        'Clipping','off','Units','inches', 'FontSize',8);
    ah.Position(3:4) = [2 2];
    
    vx = repmat([.5;.5;.25;.25;.5],1,numel(rval_labels{i}));
    vy = repmat([0;step_sz,;step_sz;0;0],1,numel(rval_labels{i}));
    vy = vy + repmat(1:step_sz:size(d,2),5,1) + repmat([.5;-.5;-.5;.5;.5],1,numel(rval_labels{i}));
    patch('XData',vx,'YData',vy,'FaceColor','k');
    title(assay_labels{i});
    axis('equal','tight');
end

%% generate an r-value plot for each assay, with a line for each metric as a function of days between

fh = gobjects(numel(rval_labels),1);
ahs = cellfun(@(r) gobjects(numel(r),1), rval_labels, 'UniformOutput', false);
for i=1:numel(rval_labels)
    fh(i) = figure('Name',assay_labels{i},'Units','inches');
    hold on;
    p = metric_permutation{i};
    d=data_by_metric{i}(p);
    for j=1:numel(rval_labels{i})
        ahs{i}(j) = subplot(2,5,j);
        nreps = 100;
        r = rval_by_metric{i}{p(j)};
        ci95 = NaN(nreps,size(r,1)-1);
        for rep=1:nreps
           bs_r = corrcoef(d{j}(randi(size(d{j},1),[size(d{j},1) 1]),:),'rows','pairwise');
           ci95(rep,:) = arrayfun(@(num_r) nanmean(diag(bs_r,num_r)), 1:size(bs_r,1)-1);
        end
        rval_line = mean(ci95);
        ci95 = [prctile(ci95,97.5);prctile(ci95,2.5)];
        xx = linspace(-1,1,numel(rval_line));
        hold on;
        vx = [xx fliplr(xx)];
        vy = [ci95(2,:) fliplr(ci95(1,:))];
        patch('XData',vx(:),'YData',vy(:),'FaceColor',[.75 .75 .75],'EdgeColor','none');
        plot(xx,rval_line,'k','LineWidth',1);
        plot([-1 1],[0 0],'k--');
        axis(ahs{i}(j),'equal');
        set(ahs{i}(j),'YLim',[-1 1],'XLim',[-1 1],'XTick',xx(2:2:end),...
            'XTickLabel',2:2:numel(rval_line),'FontSize',6,'Units','inches');
        
        xlabel('days','FontSize',7);
        ylabel('r-value','FontSize',7);
        title(pretty_labels(rval_labels{i}(p(j))),'FontSize',7);
        drawnow;
        
    end
end

step = 1.1;
for i=1:numel(rval_labels)
    for j=1:numel(rval_labels{i})
        pos = [mod(j-1,5) + 1, 3 - ceil(j/5), .75 .75];
        if (3 - ceil(j/5)) > 1
            pos(2) = pos(2) + .25;
        end
        pos(1) = pos(1) * step;
        pos(1:2) = pos(1:2)-.5;
        ahs{i}(j).Position = pos;
    end
end

%% generate an r-value plot for each assay, with a line for each metric as a function of days between

figure;
ah=subplot(2,3,1);
ah2 = subplot(2,3,2);
rp = [10,8,2,7,5,1,3,4,6,9];
d = cat(1,data_by_metric{2}(rp));
d = cat(2,d{:});
plotCorr(d,'Cluster',false,'Patch',false,'Parent',[ah ah2]);
delete(ah2);
title('LED Y-maze Correlation Matrix');
axis(ah,'equal','tight');

for i=1:numel(rval_labels)
     ah = subplot(2,3,i+1);
     hold(ah,'on');
     lhs = gobjects(numel(rval_labels{i}),1);
     for j=1:numel(rval_labels{i})
         r = rval_by_metric{i}{j};
         rval_line = NaN(size(r,1)-1,1);
         for k=1:numel(rval_line)
             rval_line(k) = nanmean(diag(r,k));
         end
         xx = linspace(-1,1,numel(rval_line));
         lhs(j) = plot(xx,rval_line,'LineWidth',1);
     end
     plot([-1 1],[0 0],'k--');
     legend(lhs,pretty_labels(rval_labels{i}),'Location','SouthWest','FontSize',6);
     axis(ah,'equal');
     set(ah,'YLim',[-1 1],'XLim',[-1 1],'XTick',xx,'XTickLabel',1:numel(rval_line));
     ylabel('r-value');
     xlabel('days between measurements');
     title(unique_assays{i});
end

