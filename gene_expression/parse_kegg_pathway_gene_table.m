
% read in raw text
p = ['C:\Users\winsl0w\Documents\decathlon\decathlon_paper_code\'...
    'decathlon\gene_expression\kegg_all_pathway_genes.txt'];

fid = fopen(p);
raw = textscan(fid,'%s','delimiter','\t');
raw = raw{1};
fclose(fid);

%%

% parse pathways and genes
all_pathways = raw(logical(mod(1:numel(raw),2)));
all_pathways = regexp(all_pathways,'(?<=:).*','match');
all_pathways = cellfun(@(s) s{1}, all_pathways, 'UniformOutput', false);
kegg_all_genes = raw(~mod(1:numel(raw),2));
kegg_all_genes = regexp(kegg_all_genes,'(?<=:).*','match');
kegg_all_genes = cellfun(@(s) s{1}, kegg_all_genes, 'UniformOutput', false);

% get unique entries
[u_pathways,~,pathway_idx] = unique(all_pathways,'stable');
n_pathways = numel(u_pathways);
u_genes = unique(kegg_all_genes,'stable');
n_genes = numel(u_genes);

% compute overlap between decathlon and kegg genes
%[fbgn,has_match] = fbgn2kegg(u_genes,'kegg');
has_match = cellfun(@(s) strcmp(D_rnaseq_models.kegg,s), u_genes, 'UniformOutput', false);
has_match = cat(2,has_match{:});
n_match = sum(any(has_match));
fprintf('%i of %i (%0.2f) modeled decathlon genes have an associated pathway\n',...
    n_match,numel(D_rnaseq_models.fbgn),n_match/numel(D_rnaseq_models.fbgn));

% parse genes by pathway
pathway_genes = cell(n_pathways,1);
for i=1:numel(u_pathways)
    pathway_genes{i} = kegg_all_genes(pathway_idx==i);
end

% calculate fraction of total genes in the decathlon set
frac_genes = NaN(n_pathways,1);
num_genes = NaN(n_pathways,1);
for i=1:numel(pathway_genes)
    num_genes(i) = sum(cellfun(@(s) any(strcmp(D_rnaseq_models.kegg,s)), pathway_genes{i}));
    frac_genes(i) = num_genes(i)/numel(pathway_genes{i});
end

% remove pathways with < 2 genes that cannot be enriched
num_filt = num_genes > 1;
num_genes = num_genes(num_filt);
frac_genes = frac_genes(num_filt);
u_pathways = u_pathways(num_filt);
pathway_genes = pathway_genes(num_filt);

% filter out pathways not hit in decathlon
[cat_names,cat_p] = unique(cat(1,D_enrichment_results.cat_labels),'stable');
dec_pathways = unique(cat(1,D_enrichment_results.cat_id),'stable');
pathway_filt = cellfun(@(s) any(strcmp(dec_pathways,s)), u_pathways);
pathway_genes = pathway_genes(pathway_filt);
pathways = u_pathways(pathway_filt);
num_genes = num_genes(pathway_filt);
frac_genes = frac_genes(pathway_filt);
fprintf('%i of %i (%0.2f) pathways hit\n',...
    sum(pathway_filt),numel(pathway_filt),sum(pathway_filt)/numel(pathway_filt));

% match cat sorting order
[dec_pathways,dec_p] = sort(dec_pathways);
cat_names = cat_names(dec_p);
[pathways,kegg_p] = sort(pathways);
pathway_genes = pathway_genes(kegg_p);
num_genes = num_genes(kegg_p);
frac_genes = frac_genes(kegg_p);

subplot(2,1,1);
histogram(frac_genes,linspace(0,1,10));
xlabel('fraction of total category genes in dataset');
set(gca,'TickLength',[0 0],'YLim',[0 20]);


% parse the decathlon genes by category
nhit_per_cat = NaN(numel(num_genes),1);
for i=1:numel(cat_names)
    g = [];
    for j=1:numel(D_enrichment_results)
        cat_idx = find(strcmp(D_enrichment_results(j).cat_labels,cat_names{i}));
        if cat_idx
            tmp_g = D_enrichment_results(j).gene_fb_gene_num(D_enrichment_results(j).prob_gene_given_cat(cat_idx,:)>0);
            tmp_g = tmp_g(~cellfun(@isempty,tmp_g));
            g = [g;tmp_g];
        end
    end
    nhit_per_cat(i) = numel(unique(g));
end

% fraction hit per cat
frac_hit_per_cat = nhit_per_cat./num_genes;
subplot(2,1,2);
histogram(frac_hit_per_cat,linspace(0,1,10));
xlabel('fraction of overlapping genes hit per category');
set(gca,'TickLength',[0 0],'YLim',[0 20]);

[ff,perm] = sort(frac_hit_per_cat,'descend');


% append gene lists to enrichment data struct
for i=1:numel(D_enrichment_results)
   cats = D_enrichment_results(i).cat_id;
   D_enrichment_results(i).cat_genes = cell(size(cats));
   for j=1:numel(cats)
      cat_idx = find(strcmp(pathways,cats{j}));
      D_enrichment_results(i).cat_genes(j) = pathway_genes(cat_idx);
   end
end


%% plot category avg. min bootstrapped p-value against pathway metrics


pvals = cat(1,D_enrichment_results.avg_min_pval);
pvals = pvals(cat_p);
pvals = pvals(dec_p);

figure;
subplot(3,1,1);
pretty_scatter(pvals,num_genes,'k');
xlabel('avg. -log[p]');
ylabel('num. overlapping genes');
title(sprintf('r=%0.2f',corr(pvals,num_genes)));
set(gca,'XLim',[0 15]);

subplot(3,1,2);
pretty_scatter(pvals,frac_genes,'k');
xlabel('avg. -log[p]');
ylabel('fraction overlapping genes');
title(sprintf('r=%0.2f',corr(pvals,frac_genes)));
set(gca,'XLim',[0 15],'YLim',[0 1]);


subplot(3,1,3);
pretty_scatter(pvals,frac_hit_per_cat,'k');
xlabel('avg. -log[p]');
ylabel('fraction overlapping genes hit');
title(sprintf('r=%0.2f',corr(pvals,frac_hit_per_cat)));
set(gca,'XLim',[0 15],'YLim',[0 1]);









