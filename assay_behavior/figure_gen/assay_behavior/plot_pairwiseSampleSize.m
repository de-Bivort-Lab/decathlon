function fig_handle = plot_pairwiseSampleSize(D)
% generate a matrix of sample sizes for pairwise comparisons of behavioral
% metrics in decathlon data struct D

% calculate number of samples for each pairwise comparison
nFields = numel(D.fields);
n = NaN(nFields);
for i=1:nFields
    for j=1:nFields       
        if i~=j
            n(i,j) = sum(~isnan(D.data(:,i)) & ~isnan(D.data(:,j)));
        end
    end
end

fig_handle = figure();
imagesc(n);
colormap('hot');
colorbar

% format field labels for display
fLabels = D.fields;
for i = 1:length(fLabels)
    tmp = fLabels{i};
    tmp(tmp=='_')=' ';
    fLabels(i)={tmp};
end

fsz = 10;
set(gca,'fontsize', fsz);
set(gca,'Ytick',1:nFields,'YtickLabel', fLabels,'fontsize',8, 'TickLength', [0 0]);
set(gca,'XTick',1:nFields,'XTickLabel', fLabels,'fontsize',8,'XTickLabelRotation',45);
title('Pairwise Sample Size');