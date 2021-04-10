
% load behavioral data
D_b = load_decathlon_structs(pwd,'D_us');
D_b(2).ID = cat(1,D_b(1).ID(end),D_b(2).ID);
D_b(2).pdfs = cat(1,D_b(1).pdfs(end,:),D_b(2).pdfs);
D_b(1).ID(end) = [];
D_b(1).pdfs(end,:) = [];

% set rnaseq preprocessing params and pre-process
min_tot_reads = 1E6;
min_rpm = 16;
molaspass=interp1([1 51 102 153 204 256],...
        [0 0 0; 0 0 .75; .5 0 .8; 1 .1 0; 1 .9 0; 1 1 1],1:256);
    
% load rnaseq data and drop decathlon-1 data
D_seq = load_decathlon_structs(pwd,'D_seq');
D_p = preprocess_rnaseq_data(D_seq, min_tot_reads, min_rpm);
D_p = D_p(2:3);

% update behavioral data struct ID no.
for i=1:numel(D_p)
    f = ismember(D_b(i).ID,D_p(i).ID);
    D_b(i).pdfs = D_b(i).pdfs(f,:);
    D_b(i).ID = D_b(i).ID(f);
    
    f = ismember(D_p(i).ID,D_b(i).ID);
    D_p(i).data = D_p(i).data(f,:);
    D_p(i).ID = D_p(i).ID(f);
end

%%

% compute p-vals of simple gene x behavior linear models
model_p_value = cell(numel(D_p),1);
r_vals = cell(numel(D_p),1);
for i=1:numel(D_p)
    model_p_value{i} = NaN(size(D_p(i).data,2),size(D_b(i).pdfs,2));
    for j=1:size(D_b(i).pdfs,2)
        fprintf('D%i, Behavior %i of %i\n',i,j,size(D_b(i).pdfs,2));
        [r_vals{i}(:,j),model_p_value{i}(:,j)] = corr(D_b(i).pdfs(:,j),D_p(i).data);
    end
end

%%

% plot clustered heatmaps of model pval matrices
for i=1:numel(model_p_value)
    figure;
    %subplot(numel(model_p_value),2,(i-1)*2+1);
    p = -log10(model_p_value{i})';
    p(isnan(p)) = 0;
    
    [~, grp_idx] = group_unsupervised_clusters;
    grp_idx = cat(2,grp_idx{:});
    [~,g_perm] = sort(nanmean(p),'descend');

    subplot(2,2,1);
    bar(sum(p(grp_idx,g_perm)>2),'EdgeColor','none','FaceColor','k');
    
    subplot(2,2,4);
    barh(sum(p(grp_idx,g_perm)>2,2),'EdgeColor','none','FaceColor','k');
%     
    subplot(2,2,3);
    imagesc(p(grp_idx,g_perm));
    colormap(flip(bone,1));
    caxis([0 3]);
    colorbar;
    title(sprintf('decathlon-%i - mean expression sort',i));
    xlabel('genes');
    ylabel('behaviors');
    set(gca,'TickLength',[0 0], 'Clipping', 'off');
    
    hold on;
    [~, grp_idx] = group_unsupervised_clusters;
    vx = [-1;-1;-.5;-.5;-1] .* 100;
    y = cellfun(@numel, grp_idx);
    y = [cumsum(y)-y+1 cumsum(y)] + [-.25 .25];
    vy = y(:,[1 2 2 1 1])';
    patch('XData',repmat(vx,1,size(vy,2)),'YData',vy,'FaceColor','r','EdgeColor','none');
end

%%

% set significance threshold and export path
fdir = uigetdir();
alpha = 0.01;
save_dir = cat(2,fdir,'/RNAseq_unsupervised_bootstrap/');
if ~exist(save_dir,'dir')
    d1 = arrayfun(@(i) sprintf('d%i',i), 1:numel(D_p), 'UniformOutput', false);
    d2 = {'shuffled';'unshuffled'};
    for i=1:numel(d1)
        for j=1:numel(d2)
            mkdir(fullfile(save_dir,'bs_data',d1{i},d2{j}));
        end
    end
    mkdir(fullfile(save_dir,'obs_data'));
    mkdir(fullfile(save_dir,'output'));
end

% lookup FBGN to KEGG gene id conversion
[bg_kegg_ids,has_kegg_gid] = convert_gene_ids(D_p(1).geneID,'fbgn','kegg');

% Export observed model hits
for i=1:numel(model_p_value)
    fname = sprintf('%sobs_data/D%i_pval_obs.csv',save_dir,i);
    [grp_name, grp_idx] = group_unsupervised_clusters;
    pvals = model_p_value{i}(has_kegg_gid,cat(2,grp_idx{:}));
    write_csv_mat(fname,pvals<alpha,grp_idx,grp_name,bg_kegg_ids);
end


% Export bootstrapped observed and shuffled model matrices for 
% ClusterProfiler analysis in R

% query field names and group aprior fields
output_dir = sprintf('%sbs_data',save_dir);

% iterate over decathlon batches and aprior groups
nreps = 1000;
for i=1:numel(D_b)
    for rep = 1:nreps
        
        fprintf('%i of %i\n',rep,nreps);
        
        idx = randi(size(D_b(i).pdfs,1),[size(D_b(i).pdfs,1) 1]);
        
        % bootstrap resample individual behavior and gene data
        [~, grp_idx] = group_unsupervised_clusters;
        b = D_b(i).pdfs(idx,cat(2,grp_idx{:}));
        g = D_p(i).data(idx,has_kegg_gid);

        % compute model pvals
        pvals = NaN(size(g,2),size(b,2));
        for j=1:size(b,2)
            [~,pvals(:,j)] = corr(b(:,j),g);
        end

        % export hits to file
        pvals = pvals < alpha;
        fname = sprintf('%s/d%i/unshuffled/D%i_unshuffled_pval_bs_%i.csv',output_dir,i,i,rep);
        write_csv_mat(fname,pvals,grp_idx,grp_name,bg_kegg_ids);
        
        % shuffle individual behavior and gene data
        b = shuffle_columns(D_b(i).pdfs(:,cat(2,grp_idx{:})));
        g = shuffle_columns(D_p(i).data(:,has_kegg_gid));

        % compute shuffled model pvals
        pvals = NaN(size(g,2),size(b,2));
        for j=1:size(b,2)
            [~,pvals(:,j)] = corr(b(:,j),g);
        end

        % export hits to file
        pvals = pvals < alpha;
        fname = sprintf('%s/d%i/shuffled/D%i_shuffled_pval_bs_%i.csv',output_dir,i,i,rep);
        write_csv_mat(fname,pvals,grp_idx,grp_name,bg_kegg_ids);
    end
end




%% save model data and meta data to struct

D_rnaseq_models_us = struct('fbgn',[],'kegg',[],'gene_names',[],...
    'gene_symbols',[],'pvals',[],'metric_labels',[]);

D_rnaseq_models_us.fbgn = D_p(1).geneID;
D_rnaseq_models_us.kegg = cell(numel(D_p(1).geneID),1);
[g,f] = convert_gene_ids(D_p(1).geneID,'fbgn','kegg');
D_rnaseq_models_us.kegg(f) = g;
D_rnaseq_models_us.gene_names = cell(numel(D_p(1).geneID),1);
[g,f] = convert_gene_ids(D_p(1).geneID,'fbgn','name');
D_rnaseq_models_us.gene_names(f) = g;
D_rnaseq_models_us.gene_symbols = cell(numel(D_p(1).geneID),1);
[g,f] = convert_gene_ids(D_p(1).geneID,'fbgn','symbol');
D_rnaseq_models_us.gene_symbols(f) = g;
D_rnaseq_models_us.pvals = model_p_value;

D_rnaseq_models_us.metric_labels = cell(size(model_p_value));
[~, grp_idx] = group_unsupervised_clusters;
for i=1:numel(model_p_value)
    grp = find(cellfun(@(j) any(ismember(i,j)), grp_idx));
    j = find(i==grp_idx{grp});
    D_rnaseq_models_us.metric_labels{i} = sprintf('%s_%i',grp_name{grp},j);
end

