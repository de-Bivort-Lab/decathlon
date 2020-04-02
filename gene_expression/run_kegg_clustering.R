# import cluster profiler
library('clusterProfiler')


# define paths
input_path = "C:/Users/winsl0w/Documents/decathlon/decathlon_paper_data/decathlon_go_geneid_data/apriori_grp_hits/kegg_gid/";
output_path = "C:/Users/winsl0w/Documents/decathlon/decathlon_paper_data/decathlon_go_geneid_data/apriori_grp_hits/kegg_enrichments/";
all_paths = list.files(input_path);

# import gene background
geneBackground <- read.csv(paste(input_path,"background_kegg_gid.csv",sep=""),sep="\t", stringsAsFactors = FALSE);
geneBackground = geneBackground[[1]];

# initialize placeholders
feature = c();
pathway = c();
p_value = c();
num_genes = c();

# iterate over gene lists
for(p in all_paths){
  
  # read in file
  geneList <- read.csv(paste(input_path,p,sep=""),sep="\t", stringsAsFactors = FALSE);
  geneList = geneList[[1]];
  
  # do kegg enrichment analysis
  kk <- enrichKEGG(gene = geneList,organism = 'dme', pvalueCutoff = 0.05);
  
  # look for hits
  results = attr(kk,"result");
  is_hit = results[[6]]<0.05;
  num_hits = sum(is_hit, na.rm = TRUE);
 
  
  if(num_hits){
    gene_hits = results[[8]][is_hit];
    p_value = append(p_value,results[[6]][is_hit])
    num_genes = append(num_genes,unlist(lapply(gene_hits, function(x) length(unlist(strsplit(x,"/"))))))
    pathway = append(pathway,results[[2]][is_hit])
    tmp_feature = strsplit(p,"\\.");
    feature = append(feature,rep(tmp_feature[[1]][1], times = num_hits))
  }
}

# construct summary data frame
kegg_summary = data.frame(feature,pathway,p_value,num_genes)

# export to file
write.csv(kegg_summary, past(output_path,"kegg_summary_metrics.csv",sep=""))






