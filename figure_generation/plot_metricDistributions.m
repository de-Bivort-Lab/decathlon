function plot_metricDistributions(D, varargin)
% plot metric kde for decathlon data struct or array of decathlon data structs D
% If D is an array of decathlon data structs, metrics will be paired
% across data structs and plotted on the same axis

% parse inputs
labels = '';
for i=1:numel(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'Labels'
                i=i+1;
                labels = varargin{i};
        end
    end
end

% ensure field names have standard format
for i=1:numel(D)
   D(i).fields = standardize_fieldnames(D(i).fields); 
end

% iterate over data structures
[nFields, idx] = max(arrayfun(@(d) numel(d.fields), D));
[D_trim,~] = arrayfun(@(d) pairFields(d,d,'Trim',true), D, 'UniformOutput', false);
D_trim = cat(1,D_trim{:});
field_list = unique(cat(1,D_trim.fields));

% iterate over fields
nRows = 8;
nCols = 4;
colors = {'b',[1 0.5 0],'m'};
plot_opts = {'Color', 'b', 'LineWidth', 0.5};
for i=1:numel(field_list)
    
    % open new figure if necessary
    if mod(i-1,nRows*nCols)==0
       figure; 
    end
    subplot(nRows,nCols,mod(i-1,nRows*nCols)+1);
    
    % get current field
    f = field_list{i};
    
    % iterate over data structs
    all_dat = [];
    all_bins = [];
    for j=1:numel(D)
        % find matching data field
        fidx = find(strcmpi(D_trim(j).fields,f),1);
        fdat = D_trim(j).data(:,fidx);
        if any(~isnan(fdat))
            plot_opts(2) = colors(j);
            [~, bins] = ...
                autoPlotDist(fdat,~isnan(fdat),gca,'PlotOptions',plot_opts);
            all_dat = [all_dat; fdat(~isnan(fdat))];
            all_bins = [all_bins bins];
            set(gca,'XLim',[min(all_bins) max(all_bins)],...
                'XTick',linspace(round(min(all_bins),1),round(max(all_bins),1),3));
            ylabel('p');
            x_label = f;
            x_label(x_label=='_') = ' ';
            xlabel(x_label);
        end
        set(gca,'XLim',[-4 4]);
    end
end
