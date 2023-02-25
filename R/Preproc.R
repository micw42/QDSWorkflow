#' @export
BulkPreproc = function(df, cell_thresh=2, count_thresh=0.5) {
  if (all(row.names(df)==seq(1, nrow(df)))) {
    gene_col = colnames(df)[1]
    df = df %>% distinct(!!sym(gene_col), .keep_all=T) %>%
      column_to_rownames(var=gene_col)
  }
  df = df %>% bulk_preproc(cell_thresh=cell_thresh, count_thresh=count_thresh)
  return(df)
}

#' @export
scPreproc = function(df, prop=0.002, 
                     nUMI_filt=500, nGene_filt=250, log10_filt=0.75, 
                     mito_filt=0.2, variable.features=10000) {
  if (all(rownames(df)==seq(1, nrow(df)))) {
    gene_col = colnames(df)[1]
    df = df %>% distinct(!!sym(gene_col), .keep_all=T) %>%
      column_to_rownames(var=gene_col)
  }
  df = df %>% sc_preproc(prop=prop, 
                         nUMI_filt=nUMI_filt, nGene_filt=nGene_filt, 
                         log10_filt=log10_filt, mito_filt=mito_filt, 
                         variable.features=variable.features)
  return(df)
}
