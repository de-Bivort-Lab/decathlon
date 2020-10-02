% Plot gene expression and gene x behavior model heatmaps

fdir = uigetdir(pwd,'Select the decathlon_paper_data directory');
D_b = load_decathlon_structs(fdir,'D_als_filled');
D_seq = load_decathlon_structs(fdir,'D_seq');
D_p = pair_rnaseq(D_seq);
pairs = unique_idx_pairs(numel(D_p),1);
min_reads = 1E5;

% --- GENE EXPRESSION HEATMAPS ----- %
figure;

% get expression data and filter out empty reads
d_norm = cell(numel(D_p),1);
for i=1:numel(D_p)
    d = D_p(i).data;
    filt = nansum(d,2)<min_reads;
    d = d(~filt,:);
    
    % normalize to reads per million
    scale_rpm = sum(d,2)./min(sum(d,2));
    d_norm{i} = d./repmat(scale_rpm,1,size(d,2));
    
    % create read count plot
    subplot(numel(D_p),1,i);
    molaspass=interp1([1 51 102 153 204 256],...
        [0 0 0; 0 0 .75; .5 0 .8; 1 .1 0; 1 .9 0; 1 1 1],1:256);
    [~,col_perm] = sort(nanmean(d_norm{i}));
    [~,row_perm] = sort(nanmean(d_norm{i},2));
    imagesc(log10(d_norm{i}(row_perm,col_perm)));
    colormap(molaspass);
    cb = colorbar;
    caxis([0 2]);
    cb.Label.String = 'log reads';
    xlabel('genes');
    ylabel('individual flies');
    title(sprintf('decathlon-%i',i));
end

% insert missing individual data into D1 struct
f = true(192,1);
f([99 114])=false;
d = NaN(size(D_b(1).data));
d(f,:) = D_b(1).data;
D_b(1).data = d;

% set rnaseq preprocessing params and pre-process
min_tot_reads = 1E6;
min_rpm = 10;
molaspass=interp1([1 51 102 153 204 256],...
        [0 0 0; 0 0 .75; .5 0 .8; 1 .1 0; 1 .9 0; 1 1 1],1:256);
    
% preprocess rnaseq data
D_p = preprocess_rnaseq_data(D_seq, min_tot_reads, min_rpm);

% ---- EXPRESSION SCREE PLOTS ---- %
figure;
for i=1:numel(D_p)
    ah = subplot(numel(D_p),1,i);
    fprintf('Bootstrapping PCA for dataset %i of %i\n',i,numel(D_p));
    plot_pca_bootstrap(D_p(i).data,100);
    set(ah,'XScale','log','XLim',[1 100]);
end

% update behavioral data struct ID no.
for i=1:numel(D_p)   
    D_b(i).data = D_b(i).data(D_p(i).ID,:);
end

% compute p-vals of simple gene x behavior linear models
model_p_value = cell(numel(D_p),1);
for i=1:numel(D_p) 
    model_p_value{i} = NaN(size(D_p(i).data,2),size(D_b(i).data,2));
    for j=1:size(D_b(i).data,2)
        fprintf('D%i, Behavior %i of %i\n',i,j,size(D_b(i).data,2));
        [~,model_p_value{i}(:,j)] = corr(D_b(i).data(:,j),D_p(i).data);
    end
end

% ---- CLUSTERED MODEL P-VAL HEATMAPS ---- %
for i=1:numel(model_p_value)
    figure;
    p = -log10(model_p_value{i})';
    p(isnan(p)) = 0;
    [~, ~, grp_idx] = group_apriori_fields(D_b(i));
    grp_idx = cat(1,grp_idx{:});
    [~,g_perm] = sort(nanmean(p),'descend');

    subplot(2,2,1);
    bar(sum(p(grp_idx,g_perm)>2),'EdgeColor','none','FaceColor','k');
    
    subplot(2,2,4);
    barh(sum(p(grp_idx,g_perm)>2,2),'EdgeColor','none','FaceColor','k');
    
    subplot(2,2,3);
    imagesc(p(grp_idx,g_perm));
    colormap(flip(bone,1));
    caxis([0 3]);
    colorbar;
    title(sprintf('decathlon-%i - mean expression sort',i));
    xlabel('genes');
    ylabel('behaviors');
    set(gca,'TickLength',[0 0]);
end


