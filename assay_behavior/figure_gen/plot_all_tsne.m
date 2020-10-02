function plot_all_tsne(data,varargin)

% set defaults
color = [];
labels = [];
plot_opt = {};
distance_type = 'euclidean';

% parse optional inputs
for i=1:numel(varargin)
    arg = varargin{i};
    if ischar(arg)
       switch arg
           case 'Color'
               i = i + 1;
               color = varargin{i};
           case 'Labels'
               i = i + 1;
               labels = varargin{i};
           case 'PlotOptions'
               i = i + 1;
               plot_opt = varargin{i};
           case 'Distance'
               i = i + 1;
               distance_type = varargin{i};
       end   
    end
end


fh = figure;
dist_type = 'euclidean';
perplex = 8;
batch = {'D12';'D3'};
leg_loc = 'NorthEast';
fdir = 'D:\decathlon_data_and_analysis\decathlon_analysis\figures\multidimensional_scaling\';
fname = sprintf('D12_D3_tsne_perplex%i_%s.fig',perplex,dist_type);
ncol = 3;
opts = statset;
opts.MaxIter = 5000;
opts.TolFun = 1E-12;

% Initialize color variables
hsv_color = ones(numel(apriori_data),3);
hue = linspace(0,1,numel(apriori_data)+1);
hsv_color(:,1) = hue(1:numel(apriori_data));
hsv_color(:,2) = rand(numel(apriori_data),1).*.75 + .25;
rgb_color = hsv2rgb(hsv_color);

D=load_decathlon_structs;
D_p = pair_decathlon_structs(D);
D_p = cat_decathlon_structs(D_p);
for i=1:numel(D_p)

    ah = subplot(2,ncol,(i-1)*ncol+1);

    % caluclate fraction missing
    pct_miss = sum(D_p(i).imputed,2)./numel(D_p(i).fields);
    pct_miss = zeros(size(D_p(i).data,1),1);
    pct_miss(383:end) = 1;
    embedded_data = tsne(D_p(i).data,'Distance',dist_type,'Perplexity',20,'Algorithm','exact','Options',opts);
    scatter(embedded_data(:,1),embedded_data(:,2),12,pct_miss,'filled');
    title(sprintf('t-SNE %s - individuals (%s)',dist_type,batch{i}));
    cb = colorbar;
    cb.Label.String = 'fraction missing';
    caxis([0 1]);

    [apriori_data, apriori_names, p] = group_apriori_fields(D_p(i));

    % concatenate grouped data
    cat_data = cat(2,apriori_data{:});
    embedded_data = tsne(cat_data','Distance',dist_type,'Perplexity',perplex,'Algorithm','exact','Options',opts);

    % scatter points and color by apriori group
    grp_idx = cellfun(@(ad,ii) repmat(ii,size(ad,2),1), apriori_data,...
        num2cell(1:numel(apriori_data))','UniformOutput', false);
    grp_idx = cat(1,grp_idx{:});

    ah = subplot(2,ncol,(i-1)*ncol+2);
    hold on;
    lh = gobjects(numel(apriori_data),1);
    for j=1:numel(apriori_data)
        grp_data = embedded_data(grp_idx == j, :);
        lh(j) = pretty_scatter(grp_data(:,1),grp_data(:,2),rgb_color(j,:),'MarkerSize',6);
    end
    title(sprintf('t-SNE %s - metrics (%s)',dist_type,batch{i}));
    legend(lh,apriori_names,'Location',leg_loc,'FontSize',6);
    buff = max([diff(ah.XLim) diff(ah.YLim)])*.2;
    set(ah,'XLim',ah.XLim + [-buff/2 buff*2],'YLim',ah.YLim + [-buff/2 buff*2]);


    % group data by assay
    [assay,~] = parse_fieldnames(D_p(i).fields(cat(1,p{:})));
    assays = unique(assay);
    grp_idx = cellfun(@(a) strcmpi(assay,a), assays, 'UniformOutput', false);
    grp_idx = sum(cat(2,grp_idx{:}).*repmat(1:numel(assays),numel(assay),1),2);

    ah = subplot(2,ncol,(i-1)*ncol+3);
    hold on;
    lh = gobjects(numel(assays),1);
    for j=1:numel(assays)
        grp_data = embedded_data(grp_idx == j, :);
        lh(j) = pretty_scatter(grp_data(:,1),grp_data(:,2),rgb_color(j,:),'MarkerSize',6);
    end
    title(sprintf('t-SNE %s - metrics (%s)',dist_type,batch{i}));
    legend(lh,assays,'Location',leg_loc,'FontSize',6);
    buff = max([diff(ah.XLim) diff(ah.YLim)])*.2;
    set(ah,'XLim',ah.XLim + [-buff/2 buff*2],'YLim',ah.YLim + [-buff/2 buff*2]);
end