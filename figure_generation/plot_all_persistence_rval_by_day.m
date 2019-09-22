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
rval_by_metric = cell(numel(unique_assays),1);
rval_labels = cell(numel(unique_assays),1);
for i=1:numel(unique_assays)
   [~,metric,day] = parse_fieldnames(assay_grouped_fields{i});
   unique_metrics = unique(metric);
   rval_by_metric{i} = cell(numel(unique_metrics),1);
   data_by_metric{i} = cell(numel(unique_metrics),1);
   rval_labels{i} = unique_metrics;
   metric_idx = str_list_contains(metric,unique_metrics);
   
   for j=1:numel(unique_metrics)
       metric_data = assay_grouped_data{i}(metric_idx(:,j));
       metric_day = day(metric_idx(:,j));
       [~,day_sort] = sort(metric_day);
       metric_data = metric_data(day_sort);
       data_by_metric{i}{j} = cat(2,metric_data{:});
       rval_by_metric{i}{j} = corr(data_by_metric{i}{j},'Type','spearman','rows','pairwise');
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

