%%
% 
% [D, als_data] = avg_als_impute(D, 100);
% 
% %% invert right bias
% 
% D_cat = [cat_decathlon_structs(D_als(1:2)) D_als(3)];
% 

D = load_decathlon_structs('D13_als');
for i=1:numel(D)
    [~,p] = groupFields(standardize_fieldnames(D(i).fields), 'right_bias');
    D(i).data(:,p) = D(i).data(:,p).*-1;
end

%%

id_group = cell(2,1);
id_assay = cell(2,1);

unique_assays = unique(parse_fieldnames(standardize_fieldnames(cat(1,D.fields))));
[~,groups,idx] = group_apriori_fields(D(1));
assay_colors = jet(numel(unique_assays));

collapse_fields = 'all';
collapse_mode = 'PCA';
D_p = remove_keyword_fields(D,'','KeywordGroup','none');
D_p = pair_decathlon_structs(D_p,'CollapseMode','average',...
    'CollapseFields','none','ImputeMode','none','Standardize',false);
D_p = collapseMetrics(D_p,'CollapseMode',collapse_mode,...
    'CollapseFields',collapse_fields,'PCs',2);
D_p = standardize_by_field(D_p);
permutations = cell(numel(D),1);
if ~strcmpi(collapse_fields,'all')
    for i=1:numel(D)
        [~,groups,idx] = group_apriori_fields(D_p(i));
        id_group{i} = arrayfun(@(ii,jj) repmat(ii,numel(jj{1}),1), ...
            1:numel(idx), idx', 'UniformOutput', false);
        idx = cat(1,idx{:});
        assay = parse_fieldnames(D_p(i).fields(idx));
        [~,~,id_assay{i}] = unique(assay,'stable');
        id_group{i} = cat(1,id_group{i}{:});
        permutations{i} = idx;
        D_p(i) = apply_cluster_permutation(idx,D_p(i));
    end
end


%
d_labels = {'D12';'D3'};
[fh,~,~,zp]=plotCorr(D_p(1).data,'Labels',D_p(1).fields,'Cluster',false,'Patch',false);
delete(fh);
rvals = cell(numel(D_p),1);
pvals = cell(numel(D_p),1);
for i=1:numel(D_p)
    [fh,rvals{i},pvals{i}]=plotCorr(D_p(i).data(:,zp),'Labels',D_p(i).fields(zp),...
        'Cluster',false,'Patch',false);
    %delete(fh(2));
    figure(fh(1));
    axis(gca,'equal','tight');
    title(sprintf('%s behavioral corr mat',d_labels{i}));
end

% linearize r-values between matrices and scatter
r_idx = upper_triangle_idx(size(rvals{1},1));
r_lin = cellfun(@(r) r(r_idx), rvals, 'UniformOutput', false);
figure;
pretty_scatter(r_lin{1},r_lin{2},'k');
set(gca,'XLim',[-.8 .8],'Ylim',[-.8 .8]);

figure;
PCARegressionCI([r_lin{1},r_lin{2}],gca,'XLim',[-.8 .8],'YLim',[-.8 .8]);

%% Find metric pairs with shared significant correlations

[p_idx,ii,jj] = upper_triangle_idx(size(pvals{1},1));
p_lin = cellfun(@(p) p(p_idx), pvals, 'UniformOutput', false);
p_lin = cat(2,p_lin{:});

% identify significant pairs in both matrices
is_sig = all(p_lin<0.05,2);
sig_p = p_lin(is_sig,:);
a_idx = ii(is_sig);
b_idx = jj(is_sig);

% sort by combined pval
ranks = 1:size(sig_p,1);
[~,a_sort] = sort(sig_p(:,1));
[~,b_sort] = sort(sig_p(:,2));
rank_a(a_sort) = ranks;
rank_b(b_sort) = ranks;
[~,p_sort] = sort(rank_a+rank_b);

% rank pairs
a_fields_ranked = D_p(1).fields(a_idx(p_sort));
b_fields_ranked = D_p(1).fields(b_idx(p_sort));






%% scatter select metric pairs (outliers in r-value scatter) from D12 and D3

lim = 9;
[~,min_idx] = min(abs(rvals{1}(:)+0.6871));
%[~,min_idx] = min(abs(rvals{1}(:)-0.3321));
[row1,col1] = ind2sub(size(rvals{1}),min_idx);
[~,min_idx] = min(abs(rvals{1}(:)-0.2475));
%[~,min_idx] = min(abs(rvals{1}(:)+0.06586));
[row2,col2] = ind2sub(size(rvals{1}),min_idx);
figure;
for i=1:numel(D_p)
    subplot(2,2,i);
    d = D_p(i).data;
    PCARegressionCI(d(:,[row1,col1]),gca,'XLim',[-lim lim],'YLim',[-lim lim]);
    axis('equal');
    set(gca,'XLim',[-lim lim],'YLim',[-lim lim]);
    title(d_labels{i});
    xlabel(D_p(i).fields(row1));
    ylabel(D_p(i).fields(col1));
    
    subplot(2,2,i+2);
    PCARegressionCI(d(:,[row2,col2]),gca,'XLim',[-lim lim],'YLim',[-lim lim]);
    axis('equal');
    set(gca,'XLim',[-lim lim],'YLim',[-lim lim]);
    title(d_labels{i});
    xlabel(D_p(i).fields(row2));
    ylabel(D_p(i).fields(col2));
end



%% create color patch boxes

if strcmp(collapse_mode,'PCA')
    fields = regexp(D_p(1).fields,'.*(?= \(PC[1-9])','match');
    fields = cat(1,fields{:});
    [groups,~,grp_ids] = unique(fields,'stable');
    id_group = cell(numel(D_p),1);
    id_group(:) = {grp_ids};
end

% create patch coordinates for groups
% create patches for a priori groups and assays
ph_group = gobjects(2,1);
ph_assay = gobjects(2,1);
vx_groups = cell(2,1);
vx_groups(:) = {repmat([0;0;2;2;0],1,numel(groups))};
vy_groups = cell(2,1);
vy_groups(:) = {NaN(5,numel(groups))};
color_groups = jet(numel(groups));
for i = 1:numel(vy_groups)
    for j = 1:numel(groups)
        idx = [find(id_group{i}==j,1) find(id_group{i}==j,1,'Last')];
        y = [idx(1) idx([2 2])+1 idx([1 1])]';
        y = y + 0.1.*[1;-1;-1;1;1];
        vy_groups{i}(:,j) = y;
    end
end

figure;
colormap('jet');
max_y = numel(id_assay{1});
vx_assays = cell(2,1);
vx_assays(:) = {repmat([2;2;4;4;2],1,numel(id_assay{1}))};
vy_assays = cell(2,1);
vy_assays(:) = {repmat([0;1;1;0;0],1,max_y) + repmat(1:max_y,5,1)};
for i=1:2
    subplot(1,2,i);
    hold on;
    patch('XData',vx_groups{i},'YData',max_y-vy_groups{i},...
        'FaceColor','flat','CData',linspace(0,1,numel(groups)));
    patch('XData',vx_assays{i},'YData',max_y-vy_assays{i},...
        'FaceColor','flat','CData',id_assay{i}./max(id_assay{i}));
    axis('equal','tight','off');
end


%%

[apriori_data, apriori_names, grp_idx] = group_apriori_fields(D(1));
data = cat(2,apriori_data{:});
p = cat(1,grp_idx{:});
ah1 = subplot(1,2,1);
ah2 = subplot(1,2,2);
fh1=figure;
ah3 = gca;
plotCorr(data,'Labels',D(1).fields(p),'Cluster',false,'Patch',false,'Parent',[ah2,ah3]);
close(fh1);

id_group = cellfun(@(gi,i) repmat(i,numel(gi),1),...
    grp_idx, num2cell(1:numel(grp_idx))', 'UniformOutput', false);
id_group = cat(1,id_group{:});

assays = parse_fieldnames(standardize_fieldnames(D(1).fields(cat(1,grp_idx{:}))));
unique_assays = unique(assays);
unique_assays(strcmpi(unique_assays,'circadian')) = [];
unique_assays = ['Circadian';unique_assays];
id_assay = cellfun(@(a) find(strcmpi(unique_assays,a)),assays);

% create patch coordinates for groups
% create patches for a priori groups and assays
vx_groups = repmat([0;0;2;2;0],1,numel(grp_idx));
vy_groups = NaN(5,numel(grp_idx));
color_groups = jet(numel(grp_idx));
for j = 1:numel(grp_idx)
    idx = [find(id_group==j,1) find(id_group==j,1,'Last')];
    y = [idx(1) idx([2 2])+1 idx([1 1])]';
    y = y + 0.2.*[1;-1;-1;1;1];
    vy_groups(:,j) = y;
end


colormap(ah1,'jet');
max_y = numel(id_assay);
vx_assays= repmat([2;2;4;4;2],1,numel(id_assay));
vy_assays = repmat([0;1;1;0;0],1,max_y) + repmat(1:max_y,5,1);
patch('XData',vx_groups,'YData',max_y-vy_groups,...
    'FaceColor','flat','CData',linspace(0,1,numel(grp_idx)),'Parent',ah1);
patch('XData',vx_assays,'YData',max_y-vy_assays,...
    'FaceColor','flat','CData',id_assay./max(id_assay),'Parent',ah1);
axis(ah1,'equal','tight','off');
axis(ah2,'equal','tight');

%% plot eigenvalue distributions

figure;
hold on;
colors = {[0 0 1];[1 0.5 0]};
bh = cell(numel(D),1);
for i=1:numel(D)
    d = D_p(i).data;
    [~,~,~,~,v_exp] = pca(d);
    bh{i} = bar(log(v_exp),'EdgeColor','none','FaceColor',colors{i},'FaceAlpha',0.5);
end
legend(cat(1,bh{:}),{'D12';'D3'});

%%

D_p2 = pair_decathlon_structs(D);
for i=1:numel(D_p2)
    D_p2(i).data(all(isnan(D_p2(i).data),2),:) = [];
end
all_data = cat(1,D_p2.data);
is_bk = false(size(all_data,1),1);
is_bk(1:size(D_p2(1).data,1)) = true;
embedded_data = tsne(all_data,'Perplexity',20,'Algorithm','exact');

figure; hold on;
lh1=pretty_scatter(embedded_data(is_bk,1),embedded_data(is_bk,2),[0 0 1]);
lh2=pretty_scatter(embedded_data(~is_bk,1),embedded_data(~is_bk,2),[1 0.5 0]);
axis('equal');
set(gca,'XLim',[-15 15],'YLim',[-15 15]);
title('t-SNE individuals');
legend([lh1 lh2],{'D12';'D3'});

% Initialize color variables
hsv_color = ones(numel(apriori_data),3);
hue = linspace(0,1,numel(apriori_data)+1);
hsv_color(:,1) = hue(1:numel(apriori_data));
hsv_color(:,2) = rand(numel(apriori_data),1).*.75 + .25;
rgb_color = hsv2rgb(hsv_color);

figure;
for i=1:numel(D_p2)
    D_p2(i).data(all(isnan(D_p2(i).data),2),:) = [];
    [apriori_data, apriori_names, p] = group_apriori_fields(D_p2(i));

    % concatenate grouped data
    cat_data = cat(2,apriori_data{:});
    embedded_data = tsne(cat_data','Distance',dist_type,'Perplexity',perplex,...
        'Algorithm','exact','Options',opts);

    % scatter points and color by apriori group
    grp_idx = cellfun(@(ad,ii) repmat(ii,size(ad,2),1), apriori_data,...
        num2cell(1:numel(apriori_data))','UniformOutput', false);
    grp_idx = cat(1,grp_idx{:});

    ah = subplot(1,2,i);
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
end

