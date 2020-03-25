
% load data
load(['D:\decathlon_data_and_analysis\decathlon_analysis\matrices\',...
    'decathlon_paper\nuisance unregressed\D123_olfaction_timeofday_added.mat']);

% define params
impute_mode = 'als';
k_folds = 50;

% pre-process data
D = impute_decathlon_structs(D,'ImputeMode',impute_mode);
D = standardize_by_field(D);


%%
D = D13_als;
D_labels = {'12';'3'};
pca_results = cell(4,numel(D));
for i = 1:numel(D)
    D(i).data(all(isnan(D(i).data),2),:) = [];
end


for j=1:1
    figure;
    for i=1:4
        switch i
            case 1
                collapse_mode = 'average';
                collapse_fields = 'none';
                title_str = 'circ expanded';
            case 2
                collapse_mode = 'average';
                collapse_fields = 'circadian';
                title_str = 'circ collapsed';
            case 3
                collapse_mode = 'average';
                collapse_fields = 'all';
                title_str = 'a priori collapsed';
            case 4
                collapse_mode = 'PCA';
                collapse_fields = 'all';
                title_str = 'a priori PCA';
        end
        opts = {'CollapseFields';collapse_fields;'CollapseMode';collapse_mode;...
            'PCs';2;'Trim';true};
        D_p = pair_decathlon_structs(D,opts{:});
        D_p = standardize_by_field(D_p);
        data = D_p(j).data;
        if size(data,1) > 192
            data = data(randperm(size(data,1),192),:);
        end
        

        % plot cross validation results
        subplot(4,2,(i-1)*2+1);
        [tr,te,pca_results{i,j}] = cross_validate_pca(data,'KFolds',k_folds,'TestSize',0.1);
        [~,mi] = min(mean(te,2));
        title(sprintf('D%i (%s), thresh = %i',j,title_str,mi)); 

        % scree plot
        subplot(4,2,(i-1)*2+2);
        title(sprintf('D%s, PC thresh = %i',D_labels{j},mi));
        [~,~,nkeep] = plot_pca_bootstrap(data,150,95,'noncummulative',[.7 .7 .7]);
        title(sprintf('D%s (%s), thresh = %i',D_labels{j},title_str,nkeep));
        set(gca,'YLim',[1E-2 1E2],'XLim',[1 size(data,2)]);
    end
end

%% MAKE SCREE PLOTS FOR ALL DECATHLON GROUPS OVERLAYED

D = D123_als;
pca_results = cell(4,numel(D));
colors = [0 0 1; 1 0.5 0; 1 0 1];
for i = 1:numel(D)
    D(i).data(all(isnan(D(i).data),2),:) = [];
end


for i=1:4
    switch i
        case 1
            collapse_mode = 'average';
            collapse_fields = 'none';
            title_str = 'circ expanded';
        case 2
            collapse_mode = 'average';
            collapse_fields = 'circadian';
            title_str = 'circ collapsed';
        case 3
            collapse_mode = 'average';
            collapse_fields = 'all';
            title_str = 'a priori collapsed';
        case 4
            collapse_mode = 'PCA';
            collapse_fields = 'all';
            title_str = 'a priori PCA';
    end
    subplot(2,2,i);
    opts = {'CollapseFields';collapse_fields;'CollapseMode';collapse_mode;...
            'PCs';2;'Trim';true};
    D_p = pair_decathlon_structs(D,opts{:});
    hh = cell(numel(D),1);
    for j=1:numel(D)
        %figure;

        D_col = D_p(j);
        D_col.data(all(isnan(D_col.data),2),:) = [];
        n = size(D_col.data,1);
        n(n>192)=192;
        D_col.data = D_col.data(randperm(size(D_col.data,1),n),:);
        
        % plot cross validation results
        % scree plot
        
        [~,~,nkeep,hh{j}] = plot_pca_bootstrap(D_col.data,150,95,'noncummulative',colors(j,:));
    end
    title(title_str);
    hh = cat(2,hh{:})';
    nn = {'1';'2';'3'};
    obs_labels = cellfun(@(i) sprintf('D%s observed',i), nn,'UniformOutput',false);
    shuff_labels = cellfun(@(i) sprintf('D%s shuffled',i), nn,'UniformOutput',false);
    set(cat(1,hh{:,2}),'LineStyle','--');
    legend(cat(1,hh{1:numel(D)*2}),[obs_labels;shuff_labels]);
    ph = cat(1,hh{:,3:4});
    delete(ph);
end

%% plot pca loadings for the first N PCs for each group

D = D123_als;
%D_labels = {'D12';'D3'};
D_labels = {'D1';'D2';'D3'};
npcs = 6;
for i = 1:numel(D)
    D(i).data(all(isnan(D(i).data),2),:) = [];
end
D_p = pair_decathlon_structs(D);
all_loadings = cell(numel(D),1);

for i=1:numel(D)
    figure;
    [c,s] = pca(D_p(i).data);
    all_loadings{i} = c;
    for j=1:npcs
       subplot(1,npcs,j);
       [cc,p] = sort(c(:,j));
       
       barh(1:numel(cc),cc);
%        if j==1
%            title(sprintf('PCA loadings - %s',D_labels{i}));
%        end
%        %set(gca,'XTick',[],'XLim',[0 size(c,2)+1]);
%        if j==npcs
%           set(gca,'XTick',1:size(c,2),'XTickLabels',D_p(i).fields,...
%               'XTickLabelRotation',90,'FontSize',7); 
%        end
        set(gca,'YTick',1:numel(cc),'YTickLabels',D_p(i).fields(p),'FontSize',6);
        title(sprintf('PC%i',j));
    end
end

% plot all difference matrices
pairs = unique_idx_pairs(numel(D),1);
diff_mats = cell(size(pairs,1),1);
for i=1:size(pairs,1)
    diff_mats{i} = all_loadings{pairs(i,1)} - all_loadings{pairs(i,2)};
end
clim = cellfun(@(d) [min(d(:)) max(d(:))], diff_mats, 'UniformOutput', false);
clim = cat(1,clim{:});
clim = [min(clim(:,1)) max(clim(:,2))];
figure;
for i=1:numel(diff_mats)
   subplot(1,numel(diff_mats),i);
   imagesc(diff_mats{i});
   caxis(clim);
   title(sprintf('loadings, %s - %s',D_labels{pairs(i,1)},D_labels{pairs(i,2)}));
end
figure;
d_loads = NaN(numel(diff_mats),size(c,2));
for i=1:numel(diff_mats)
   d_loads(i,:) = sqrt(sum(diff_mats{i}.^2));
end
figure;
bar(d_loads(:,1:npcs*2)');
