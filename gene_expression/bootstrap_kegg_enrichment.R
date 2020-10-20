# Run decathlon RNAseq KEGG enrichment analysis


# import cluster profiler
library('clusterProfiler')
source("decathlon/decathlon_paper_code/gene_expression/decathlon_enrich_utils.R")

# EDIT FILE PATH to set the current working directory to RNAseq_bootstrap folder
input_path = choose.dir(default = "", caption = "Select RNAseq_bootstrap folder");
setwd(input_path);

# define paths
output_path = "output/";
obs_paths = list.files("obs_data/",full.names=TRUE)
bs_dirs = list.dirs("bs_data/", full.names = TRUE, recursive = FALSE)


# initialize placeholders to store data
d12_counts = list();
d12_cats = list();
d12_gene_mats = list();
d12_pvals = list();
d12_num_metrics = list();
d12_metric_mats = list();


for(k in 1:3){
  
  # compute observed enrichment and filter duplicate results
  obs_d <- read.csv(obs_paths[k],sep=",", stringsAsFactors = FALSE);
  obs_results = concat_kegg_enrichments(obs_d)
  obs_results = filter_enrich_cats(obs_results)
  obs_cats = obs_results[[2]]
  
  # get unique list of genes
  obs_cat_genes = lapply(obs_results$geneID, function(x) unlist(strsplit(x,"/")))
  obs_all_genes = unique(unlist(obs_cat_genes))
  
  # get apriori groups
  apgs = group_apriori_fields(obs_d)
  upgs = append("all",unique(apgs[[1]]))
  
  # get subdirectory names
  bs_sub_dirs = list.dirs(bs_dirs[k], full.names = TRUE, recursive = FALSE)
  
  for(j in 1:length(bs_sub_dirs)){
    
    # get subdir description
    fdesc = strsplit(bs_sub_dirs,"/")[[j]][length(strsplit(bs_sub_dirs,"/")[[j]])]
    
    # define bootstrapped data paths
    bs_paths = list.files(bs_sub_dirs[j], full.names = TRUE);
    
    #bs_paths = bs_paths[1:5]
 
    # initialize matrices
    bs_pval_mat = replicate(length(upgs),list(matrix(1, nrow = length(obs_cats), ncol = length(bs_paths))))
    bs_num_metrics = replicate(length(upgs),list(matrix(0, nrow = length(obs_cats), ncol = length(bs_paths))))
    bs_gene_mat = replicate(length(upgs),list(matrix(0, nrow = length(obs_cats), ncol = length(obs_all_genes))))
    bs_metric_mat = replicate(length(upgs),list(matrix(0, nrow = length(obs_cats), ncol = ncol(obs_d)-1)))
    cat_counts = replicate(length(upgs),list(replicate(length(obs_cats), 0)))
    
    
    for(i in 1:length(bs_paths)){
      
      
      print(sprintf("d%i_%s, iter = %i of %i",k,fdesc,i,length(bs_paths)))
      
      # read in model pvalue data
      #bs_d <- read.csv(bs_paths[1],sep=",", stringsAsFactors = FALSE);
      bs_d <- read.csv(bs_paths[i],sep=",", stringsAsFactors = FALSE);
      
      for(L in 1:length(upgs)){
        
          # compute enrichment and restrict to observed data cats
          bs_results = concat_kegg_enrichments(bs_d, apriori_grp = upgs[L])
          bs_results = filter_enrich_cats(bs_results, filt_mode = "cat", thresh=1, cat_names = obs_cats)
        
          # get enrichments p<0.05
          p_thresh = bs_results$p.adjust < 0.05;
          metric_idx = bs_results$metric_idx[p_thresh];
          metric_cats = bs_results$Description[p_thresh];
          
          # get genelist for each enrichment category
          overlap_cats = unique(metric_cats[is.element(metric_cats,obs_cats)])
          for(cat in overlap_cats){
            # increment metric counts
            cat_mask = obs_cats == cat
            metric_mask = metric_idx[metric_cats == cat]
            bs_metric_mat[[L]][cat_mask, metric_mask] = bs_metric_mat[[L]][cat_mask, metric_mask] + 1
          }
          
          # compute avg. min p-value and avg. num metrics
          tmp_cats = bs_results$Description
          bs_cats = unique(tmp_cats)
          overlap_cats = bs_cats[is.element(bs_cats,obs_cats)]
          for(cat in overlap_cats){
            bs_pval_mat[[L]][obs_cats == cat,i] = min(bs_results$p.adjust[tmp_cats == cat])
            bs_num_metrics[[L]][obs_cats == cat,i] = sum(bs_results$p.adjust[tmp_cats == cat] < 0.05)
          }
          
          # further restrict enrichment hits to overlapping categories with p<0.05
          bs_results = filter_enrich_cats(bs_results)
          
          # count enrichment hits
          tmp_cats = bs_results$Description
          bs_cats = unique(tmp_cats)
          cat_counts[[L]] = cat_counts[[L]] + is.element(obs_cats,bs_cats)
 
          # get genelist for each enrichment category
          bs_cat_genes = lapply(bs_results$geneID, function(x) unlist(strsplit(x,"/")))
          overlap_cats = bs_cats[is.element(bs_cats,obs_cats)]
          for(cat in overlap_cats){
            overlap_cat_genes = unique(unlist(bs_cat_genes[tmp_cats == cat]))
            gene_mask = is.element(obs_all_genes,overlap_cat_genes)
            cat_mask = obs_cats == cat
            bs_gene_mat[[L]][cat_mask, gene_mask] = bs_gene_mat[[L]][cat_mask, gene_mask] + 1
          }
      }
    }
    
    # write all results to file
    for(L in 1:length(upgs)){
      
      # store counts in all counts list
      d12_counts = append(d12_counts, list(cat_counts[[L]]/length(bs_paths)))
      d12_cats = append(d12_cats, list(obs_cats))
      d12_gene_mats = append(d12_gene_mats, list(bs_gene_mat[[L]]))
      d12_num_metrics = append(d12_num_metrics, list(bs_num_metrics[[L]]))
      d12_metric_mats = append(d12_metric_mats, list(bs_metric_mat[[L]]))
      d12_pvals = append(d12_pvals, list(bs_pval_mat))
      
      dframe <- data.frame(bs_gene_mat[[L]], row.names = obs_cats)
      colnames(dframe) <- obs_all_genes
      dframe$ID = obs_results$ID
      fname = paste(output_path,"d",k,"_",fdesc,"_",upgs[L],"_cat_x_gene_mat.csv",sep="")
      write.csv(dframe, fname)
      
      dframe <- data.frame(bs_pval_mat[[L]], row.names = obs_cats)
      fname = paste(output_path,"d",k,"_",fdesc,"_",upgs[L],"_pval_mat.csv",sep="")
      write.csv(dframe, fname)
      
      dframe <- data.frame(bs_num_metrics[[L]], row.names = obs_cats)
      fname = paste(output_path,"d",k,"_",fdesc,"_",upgs[L],"_num_metrics_mat.csv",sep="")
      write.csv(dframe, fname)
      
      dframe <- data.frame(bs_metric_mat[[L]]/length(bs_paths), row.names = obs_cats)
      fname = paste(output_path,"d",k,"_",fdesc,"_",upgs[L],"_cat_x_metric_mat.csv",sep="")
      write.csv(dframe, fname)
    }
  }
}
