

D = D123_als;
max_pcs = 50;
mode = 'unshuffled';
save_fig_dir = uigetdir('Select directory to save figure outputs');

%%
for i=4
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
        'PCs';2;'Trim';true;'ImputeMode';'none'};
    D_p = D;
    D_p = pair_decathlon_structs(D_p,opts{:});
     D_p = standardize_by_field(D_p);
    
    % PCA data
    pcs = cell(numel(D_p),1);
    for j=1:numel(D_p)
        d = D_p(j).data;
        d(all(isnan(d),2),:) = [];
        d = zscore(d);
        switch mode
            case 'shuffled'
                d = shuffle_columns(d);
        end
        [pcs{j},~,~,~,v_exp] = pca(d,'NumComponents',size(d,2));
        if size(pcs{j},1) > max_pcs
            pcs{j} = pcs{j}(:,1:max_pcs);
        end
    end
    
    % iterate over unique pairs
    fh = figure('Name',sprintf('%s - %s',title_str,mode));
    pair_idx = unique_idx_pairs(numel(D_p),1);
    for j=1:size(pair_idx,1)
        % plot A dot B
        subplot(2,size(pair_idx,1),j);
        plot_pca_dot_dist(pcs{pair_idx(j,1)},pcs{pair_idx(j,2)});
        ylabel('PC-A');
        xlabel('PC-B (ranked)');
        cb = colorbar;
        ylabel(cb,'$|A \cdot B|$','Interpreter','latex');
        set(gca,'XTick',[],'YTick',[]);
        title(sprintf('A = D%i, B = D%i',pair_idx(j,1),pair_idx(j,2)));
        axis('equal','tight');
        caxis([0 0.5]);
        
        % plot B dot A
        subplot(2,size(pair_idx,1),j+size(pair_idx,1));
        plot_pca_dot_dist(pcs{pair_idx(j,2)},pcs{pair_idx(j,1)});
        ylabel('PC-A');
        xlabel('PC-B (ranked)');
        cb = colorbar;
        ylabel(cb,'$|A \cdot B|$','Interpreter','latex');
        set(gca,'XTick',[],'YTick',[]);
        title(sprintf('A = D%i, B = D%i',pair_idx(j,2),pair_idx(j,1)));
        axis('equal','tight');
        caxis([0 0.5]);
    end
    save_fig_path = sprintf('%s\\%s_%s.fig',save_fig_dir,title_str,mode);
    savefig(fh,save_fig_path);
end

%% bootstrap and average images, display delta


D = D123_als;
max_pcs = 50;
mode = 'unshuffled';
nreps = 50;

for i=4
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
        'PCs';2;'Trim';true;'ImputeMode';'none'};
    D_p = D;
    D_p = pair_decathlon_structs(D_p,opts{:});
    D_p = standardize_by_field(D_p);
    
    % PCA data
    pcs = cell(numel(D_p),nreps);
    for k=1:nreps
        fprintf('bootstrap rep %i of %i\n',k,nreps);
        for j=1:numel(D_p)
            d = D_p(j).data;
            d(all(isnan(d),2),:) = [];
            d = d(randi(size(d,1),[size(d,1) 1]),:);
            d = zscore(d);
            switch mode
                case 'shuffled'
                    d = shuffle_columns(d);
            end
            [pcs{j,k},~,~,~,v_exp] = pca(d,'NumComponents',size(d,2));
            if size(pcs{j,k},1) > max_pcs
                pcs{j,k} = pcs{j,k}(:,1:max_pcs);
            end
        end
    end
    
    % iterate over unique pairs
    pair_idx = unique_idx_pairs(numel(D_p),1);
    dp = cell(numel(D_p),nreps,2);
    for k=1:nreps
        fprintf('bootstrap rep %i of %i\n',k,nreps);
        for j=1:size(pair_idx,1)
            % plot A dot B
            dp{j,k,1} = plot_pca_dot_dist(pcs{pair_idx(j,1)},pcs{pair_idx(j,2)});

            % plot B dot A
            dp{j,k,2} = plot_pca_dot_dist(pcs{pair_idx(j,2)},pcs{pair_idx(j,1)});
        end
    end
    
    raw_avg_img = cell(size(pair_idx,1)*2,1);
    for j=1:size(pair_idx,1)
        raw_avg_img{j} = mean(cat(3,dp{j,:,1}),3);
        raw_avg_img{j+size(pair_idx,1)} = mean(cat(3,dp{j,:,2}),3);
    end
    
    % PCA data
    pcs = cell(numel(D_p),nreps);
    mode='shuffled';
    for k=1:nreps
        fprintf('bootstrap rep %i of %i\n',k,nreps);
        for j=1:numel(D_p)
            d = D_p(j).data;
            d(all(isnan(d),2),:) = [];
            d = d(randi(size(d,1),[size(d,1) 1]),:);
            d = zscore(d);
            switch mode
                case 'shuffled'
                    d = shuffle_columns(d);
            end
            [pcs{j,k},~,~,~,v_exp] = pca(d,'NumComponents',size(d,2));
            if size(pcs{j,k},1) > max_pcs
                pcs{j,k} = pcs{j,k}(:,1:max_pcs);
            end
        end
    end
    
    % iterate over unique pairs
    pair_idx = unique_idx_pairs(numel(D_p),1);
    dp = cell(numel(D_p),nreps,2);
    for k=1:nreps
        fprintf('bootstrap rep %i of %i\n',k,nreps);
        for j=1:size(pair_idx,1)
            % plot A dot B
            dp{j,k,1} = plot_pca_dot_dist(pcs{pair_idx(j,1)},pcs{pair_idx(j,2)});

            % plot B dot A
            dp{j,k,2} = plot_pca_dot_dist(pcs{pair_idx(j,2)},pcs{pair_idx(j,1)});
        end
    end
    
    shuf_avg_img = cell(size(pair_idx,1)*2,1);
    for j=1:size(pair_idx,1)
        shuf_avg_img{j} = mean(cat(3,dp{j,:,1}),3);
        shuf_avg_img{j+size(pair_idx,1)} = mean(cat(3,dp{j,:,2}),3);
    end
end

figure;
cmap=interp1([1 128 129 256],[0 0.35 0.6; .9921 1 1; 1 .9921 .9921; 0.7 0 0],1:256);
for i=1:numel(raw_avg_img)
   subplot(2,3,i);
   imagesc(raw_avg_img{i}-shuf_avg_img{i});
   colormap(cmap);
   colorbar;
   caxis([-.5 .5]);
end

figure;
for i=1:numel(raw_avg_img)
   subplot(2,3,i);
   imagesc(raw_avg_img{i});
   colorbar;
   caxis([-1 1]);
end

figure;
for i=1:numel(raw_avg_img)
   subplot(2,3,i);
   imagesc(shuf_avg_img{i});
   colorbar;
   caxis([-1 1]);
end
