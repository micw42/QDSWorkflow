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
  df = invisible(sc_preproc(df, prop=prop, 
                         nUMI_filt=nUMI_filt, nGene_filt=nGene_filt, 
                         log10_filt=log10_filt, mito_filt=mito_filt, 
                         variable.features=variable.features))
  return(df)
}

#' Filter doublets.
#' This function uses the DoubletFinder package to find doublets
#' and removes the doublets from the data table
#' 
#' @param df Data frame with raw counts. The rows should be genes, and columns should be samples. 
#' @param ann_df Metadata table containing the sample each cell belongs to. The cells are grouped by sample, and the doublet-finding algorithm is run on each group separately 
#' @param split.by Name of column in ann_df to group the cells by. It should be the column corresponding to the sample the cell was from
#' @param id_col Name of column in ann_df containing the cell IDs
#' @return Raw counts table with doublets filtered
#' @export
FilterDoublets = function(df, ann_df, split.by, byvar) {
  if (all(rownames(df)==seq(1, nrow(df)))) {
    gene_col = colnames(df)[1]
    df = df %>% distinct(!!sym(gene_col), .keep_all=T) %>%
      column_to_rownames(var=gene_col)
  }
  singlets = df %>% get_seurat_obj() %>% filter_doublets(ann_df=ann_df, split.by=split.by, byvar=byvar) %>% invisible()
  print(head(singlets))
  df = df %>% select(one_of(singlets))
  return(df)
}
