function num_active = estimate_numactive(data,fields)

[~,~,days] = parse_fieldnames(fields);
unique_days = unique(days);
num_active = NaN(numel(unique_days),1);
labels = cell(numel(unique_days),1);

for i=1:numel(unique_days)
    day = unique_days(i);
    day_idx = days==day;
    num_active(i) = sum(any(~isnan(data(:,day_idx)),2));
    labels{i}  = sprintf('Day (%i)',day);
end

bar(num_active,1);
set(gca,'XTick',1:numel(num_active),'XTickLabels',labels,...
    'XTickLabelRotation',90,'YLim',[0 size(data,1)],'TickLength',[0 0]);