function fig_handle = plot_numActive(D)

num_active = NaN(10,1);
for i=1:10
    [~,day_idx] = groupFields(D.fields,sprintf('(%i)',i));
    day_dat = D.data(:,day_idx);
    num_active(i) = sum(any(~isnan(day_dat),2));
end

% plot num active
fig_handle = figure;
plot(num_active,'LineWidth',2);
set(gca,'YLim',[0 size(D.data,1)]);
xlabel('Day');
ylabel('Number Active');
title('Flies Active');