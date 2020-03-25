function scatter_field_pairs(D, pair_idx)
% scatter pairs of fields vs each other

% calculate num row/col for subplots
n_pairs = size(pair_idx,1)*numel(D);
nCol = ceil(sqrt(n_pairs));
nCol(nCol>8) = 8;
while mod(nCol,numel(D)) ~= 0
    nCol = nCol + 1;
end
nRows = ceil(n_pairs/nCol);
nRows(nRows>8) = 8;


% initialize plot options
opts = {'Marker'; 'o'; 'LineStyle'; 'none';...
    'MarkerEdgeColor'; 'none';'MarkerSize'; 1.5; 'LineWidth'; 1};
fopts = {'FontSize',6,'FontWeight','normal'};
try
    [assay,metric] = parse_fieldnames(D(1).fields);
    metric = pretty_labels(metric);
    if numel(assay) ~= numel(metric)
        labels = D(1).fields;
    else
        labels = cellfun(@(a,m) {a;m}, assay, metric, 'UniformOutput', false);
    end
catch
    labels = D(1).fields;
end

for j=1:size(pair_idx,1)
    for i=1:numel(D)
        
        plot_num = (j-1)*numel(D) + i-1;
        subplot_num = mod(plot_num,nCol*nRows);
        if subplot_num == 0
            figure;
        end
    
        ah = subplot(nRows,nCol,subplot_num+1);

        % fit model
        x = D(i).data(:,pair_idx(j,1));
        y = D(i).data(:,pair_idx(j,2));
        
%         % scatter rank instead of value
%         [~,p] = sort(x);
%         r = 1:numel(x);
%         x(p) = r;
%         x = zscore(x);
%         [~,p] = sort(y);
%         y(p) = r;
%         y = zscore(y);
        
        x(x==0) = eps;
        y(y==0) = eps;
%         x = zscore(real(log10(x)));
%         y = zscore(real(log10(y)));
        x = zscore(real(x));
        y = zscore(real(y));
        mask = ~isnan(x) & ~isnan(y);
        x = x(mask);
        y = y(mask);

        mdl = PCARegressionCI([x y], ah);
        drawnow limitrate
        
        xlabel(labels{pair_idx(j,1)}, fopts{:});
        ylabel(labels{pair_idx(j,2)}, fopts{:});
        title(sprintf('D%i - rank %i',i,j),fopts{:});
        axis(ah,'equal');
        r=3.5;
        set(ah,'XTick',[],'YTick',[],'Clipping','on','XLim',[-r r],'YLim',[-r r]);
    end
end