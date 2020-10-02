library('stringi')

# define function to query unique apriori group names
group_apriori_fields <- function(d){
  
  # get apriori groups (exclude first column w/ KEGG gids)
  raw_names = names(d)[2:length(names(d))]
  raw_names = lapply(raw_names, function(x) unlist(strsplit(x,"[.]")))
  raw_filt = lapply(raw_names, function(x) !grepl("^[[:digit:]]",x))
  all_grp_names = list()
  for(j in 1:length(raw_names)){
    all_grp_names = append(all_grp_names, stri_paste_list(list(raw_names[[j]][raw_filt[[j]]]), sep=" "))
  }
  
  # restrict to unique list of group names and record column indices of each grp
  unique_grp_names = unlist(unique(all_grp_names))
  grp_idx = list()
  for(name in unique_grp_names){
    grp_idx = append(grp_idx, list(which(is.element(all_grp_names,name))))
  }
  
  return(list(unique_grp_names, grp_idx))
}

# filter results by lowest pvals for each unique category in enrichment list
filter_enrich_cats <- function(results, filt_mode = "pvalue", thresh = 0.05, cat_names = ""){
  
  results = results[results[[6]]<thresh,]
  cats = results[,2]
  
  if(filt_mode == "pvalue"){
    # get unique category names
    unique_cats = unique(cats);
    min_p_idx = c();
    
    # iterate over all unique categories and record index of lowest pval
    for(cat in unique_cats){
      is_cat = which(cats == cat)
      min_p_idx = append(min_p_idx,is_cat[which.min(results[is_cat,6])])
    }
    
    # output filtered category list
    return(results[min_p_idx,])
    
  } else if(filt_mode == "cat") {
    # output results matching input category list
    return(results[is.element(cats,cat_names),])
  }
}


# define function to compute kegg enrichment for specified apriori grp
concat_kegg_enrichments <- function(d, apriori_grp = "all"){

  # record background gene list
  geneBackground = d[[1]]
  
  # get apriori grps
  apgs = group_apriori_fields(d)
  if(apriori_grp == "all"){
    grp_idx = unlist(apgs[[2]]) + 1
  } else {
    grp_idx = unlist(apgs[[2]][apgs[[1]]==apriori_grp]) + 1
  }
    
  # select subset of the thresholded pval matrix corresponding to current apriori grp
  d_apriori = d[,grp_idx]
  
  # initialize placeholders
  results_list = list();
  geneset_list = list();
  
  # get number of columns in apriori grp
  if(!is.null(ncol(d_apriori))){
    nc = ncol(d_apriori)
  } else {
    nc = 1
  }
  
  # iterate over cols
  for(j in 1:nc){
    
    # compute kegg enrichment for col
    geneList = geneBackground[d_apriori[[j]]>0]
    kk <- enrichKEGG(gene = geneList,organism = 'dme', universe = geneBackground, pvalueCutoff = 0.05, keyType = "kegg");
    
    # threshold pval
    results = attr(kk,"result");
    if(length(results)){
      results$metric_idx <- replicate(nrow(results), j)
    }
    geneset = attr(kk,"geneSets");
    
    # append results to list
    results_list[[j]] = results;
    geneset_list = append(geneset_list,geneset);
  }
  
  ## concatenate results and geneset dataframes
  concat_results <- do.call("rbind", results_list);
  return(concat_results)
}
