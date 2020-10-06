% set decathlon_final_data.mat directory and load data
fdir = '';
D = load_decathlon_structs(fdir,'D13_als');

% pair structs to match fields
opts = {'CollapseFields','all','CollapseMode','PCA'};
D_p = pair_decathlon_structs(D,opts{:});

% group data by relevant fields
[grp_names,grp_idx] = groupFields(D_p(1).fields,{'clumpiness';'switchiness'});

% remove bout clumpiness measures
mask = str_list_contains(grp_names,'Bout clumpiness');
grp_names(mask) = [];
grp_idx(mask) = [];

% extract clumpiness/switchiness data and field labels
clump_switch_data = D_p(1).data(:,grp_idx);
clump_switch_labels = D_p(1).fields(grp_idx);

% define significant pairs (p < 0.005)
[~,p] = corr(clump_switch_data,'Type','spearman','rows','pairwise');
unique_p_idx = upper_triangle_idx(size(p,1));
is_sig = p(unique_p_idx)<0.005;
[r,c] = ind2sub(size(p),unique_p_idx(is_sig));
clump_switch_sig_pairs = [grp_idx(r)' grp_idx(c)'];
scatter_field_pairs(D_p(1),clump_switch_sig_pairs);

% plot loadings for those PCs
D_sig = D_p(1);
idx = unique(clump_switch_sig_pairs(:));
D_sig.data = D_sig.data(:,idx);
D_sig.fields = D_sig.fields(idx);
D_sig.loadings = D_sig.loadings(idx);
D_sig.loadings_labels = D_sig.loadings_labels(idx);
labels = D_sig.loadings_labels;
for i = 1:numel(labels)
   [assay,metric] = parse_fieldnames(labels{i});
   metric(str_list_contains(metric,'clumpiness')) = {'clump'};
   metric(str_list_contains(metric,'hand_switchiness')) = {'hand switch'};
   metric(str_list_contains(metric,'light_switchiness')) = {'light switch'};
   assay(str_list_contains(assay,'Temporal Phototaxis')) = {'T. Phototaxis'};
   assay(str_list_contains(assay,'Slow Phototaxis')) = {'S. Phototaxis'};
   labels{i} = cellfun(@(a,m) sprintf('%s %s',a,m),assay,metric,'UniformOutput', false);
end
D_sig.loadings_labels = labels;
plot_pca_loadings(D_sig);

% configure axes
ah = findall(groot,'Type','axes');
set(ah,'FontSize',6,'Units','inches');
pos = get(ah,'Position');
for i=1:numel(pos)
    pos{i}(3) = 0.5;
    pos{i}(4) = 0.75;
    ah(i).Position = pos{i};
    ah(i).Title.HorizontalAlignment = 'right';
    ah(i).Title.Position(1) = ah(i).XLim(2);
    ah(i).XLabel.String = 'feature weight';
end

% plot correlation and p-value matrices
fh = figure('Name','clump_switch corr matrix');
ah = subplot_array(2);
plotCorr(clump_switch_data,'PvalPatch',false,'Labels',clump_switch_labels,...
    'Cluster',false,'Parent',ah,'PvalPlot',true);
axis(ah,'equal','tight');






