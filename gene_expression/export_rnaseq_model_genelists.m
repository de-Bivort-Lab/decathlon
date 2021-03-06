
% load behavioral data
fdir = uigetdir(pwd,'Select the decathlon_paper_data directory');
D_b = load_decathlon_structs(fdir,'D_als_filled');

% get field sorting order
for i=1:numel(D_b)
    fields = D_b(i).fields;
    [a,m,d] = parse_fieldnames(fields);
    fields = cellfun(@(a,m) sprintf('%s %s',a,m),a,m,'UniformOutput',false);
    is_circ = strcmp(a,'Circadian');
    fields(is_circ) = cellfun(@(s,d) sprintf('%s (%i)',s,d),...
        fields(is_circ),num2cell(d(is_circ)),'UniformOutput',false);
    D_b(i).fields = fields;
    
    % sort by apriori group for d1, then match d2 to d1
    [~, grp_name, grp_idx] = group_apriori_fields(D_b(i));
    D_b(i).fields = D_b(i).fields(cat(1,grp_idx{:}));
    D_b(i).data = D_b(i).data(:,cat(1,grp_idx{:}));
end


% insert missing individual data into D1 struct
f = true(192,1);
f([99 114])=false;
d = NaN(size(D_b(1).data));
d(f,:) = D_b(1).data;
D_b(1).data = d;

% set rnaseq preprocessing params and pre-process
min_tot_reads = 1E6;
min_rpm = 16;
molaspass=interp1([1 51 102 153 204 256],...
        [0 0 0; 0 0 .75; .5 0 .8; 1 .1 0; 1 .9 0; 1 1 1],1:256);
    
% load rnaseq data
load(cat(2,fdir,'/decathlon_rnaseq_data.mat'));
D_p = preprocess_rnaseq_data(D_seq, min_tot_reads, min_rpm);

% update behavioral data struct ID no.
for i=1:numel(D_p)   
    D_b(i).data = D_b(i).data(D_p(i).ID,:);
end

% compute p-vals of simple gene x behavior linear models
model_p_value = cell(numel(D_p),1);
r_vals = cell(numel(D_p),1);
for i=1:numel(D_p) 
    model_p_value{i} = NaN(size(D_p(i).data,2),size(D_b(i).data,2));
    for j=1:size(D_b(i).data,2)
        fprintf('D%i, Behavior %i of %i\n',i,j,size(D_b(i).data,2));
        [r_vals{i}(:,j),model_p_value{i}(:,j)] = corr(D_b(i).data(:,j),D_p(i).data);
    end
end

%%

% plot clustered heatmaps of model pval matrices
for i=1:numel(model_p_value)
    figure;
    %subplot(numel(model_p_value),2,(i-1)*2+1);
    p = -log10(model_p_value{i})';
    p(isnan(p)) = 0;
    
    [~, ~, grp_idx] = group_apriori_fields(D_b(i));
    grp_idx = cat(1,grp_idx{:});
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
    set(gca,'TickLength',[0 0]);
end

% set significance threshold and export path
alpha = 0.05;
save_dir = cat(2,fdir,'/RNAseq_bootstrap/');
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
    [~, ~, grp_idx] = group_apriori_fields(D_b(i));
    pvals = model_p_value{i}(has_kegg_gid,cat(1,grp_idx{:}));
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
        
        idx = randi(size(D_b(i).data,1),[size(D_b(i).data,1) 1]);
        
        % bootstrap resample individual behavior and gene data
        [~, ~, grp_idx] = group_apriori_fields(D_b(i));
        b = D_b(i).data(idx,cat(1,grp_idx{:}));
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
        b = shuffle_columns(D_b(i).data(:,cat(1,grp_idx{:})));
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

D_rnaseq_models = struct('fbgn',[],'kegg',[],'gene_names',[],...
    'gene_symbols',[],'pvals',[],'metric_labels',[]);

D_rnaseq_models.fbgn = D_p(1).geneID;
D_rnaseq_models.kegg = cell(numel(D_p(1).geneID),1);
[g,f] = convert_gene_ids(D_p(1).geneID,'fbgn','kegg');
D_rnaseq_models.kegg(f) = g;
D_rnaseq_models.gene_names = cell(numel(D_p(1).geneID),1);
[g,f] = convert_gene_ids(D_p(1).geneID,'fbgn','name');
D_rnaseq_models.gene_names(f) = g;
D_rnaseq_models.gene_symbols = cell(numel(D_p(1).geneID),1);
[g,f] = convert_gene_ids(D_p(1).geneID,'fbgn','symbol');
D_rnaseq_models.gene_symbols(f) = g;
D_rnaseq_models.pvals = model_p_value;

D_rnaseq_models.metric_labels = cell(size(model_p_value));
for i=1:numel(model_p_value)
    D_rnaseq_models.metric_labels{i} = D_b(i).fields;
end

