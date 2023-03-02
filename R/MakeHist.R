#' Make histograms of QC metric distributions
#' 
#' 
#' The Make_Hist() function calculates total raw counts (nUMI), total number of non zero genes (nGene), 
#' non genes per total raw counts (log10GenesPerUMI), and mitochondrial gene counts per total raw counts (mitoRatio)
#'  and plots a histogram of each.
#'  
#'  
#' @param df Data frame with raw counts. The rows should be genes, and columns should be samples
#' @param ann_df Metadata table containing biological condition of each cell (optional)
#' @param grouping Name of variable to group the histograms by (optional)
#' @param mt.pattern regex pattern for mitochondrial genes (default: "^MT-")
#' 
#' @return List of histograms (ggplot2 objects) for nUMI, nGene, log10GenesPerUMI, and mitoRatio
#' @export
MakeHist = function(df, ann_df=NULL, mt.pattern="^MT-", grouping=NULL) {
  if (all(row.names(df)==seq(1, nrow(df)))) {
    gene_col = colnames(df)[1]
    df = df %>% distinct(!!sym(gene_col), .keep_all=T) %>%
      column_to_rownames(var=gene_col)
  }
  
  obj = df %>% get_seurat_obj(mt.pattern=mt.pattern)
  p1 = make_hist(obj, ann_df = ann_df, grouping=grouping, metric="nUMI")
  p2 = make_hist(obj, ann_df = ann_df, grouping=grouping, metric="nGene")
  p3 = make_hist(obj, ann_df = ann_df, grouping=grouping, metric="log10GenesPerUMI")
  p4 = make_hist(obj, ann_df = ann_df, grouping=grouping, metric="mitoRatio")
  p_list = list(p1, p2, p3, p4)
  return(p_list)
}
