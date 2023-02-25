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
