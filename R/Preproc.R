#' Filter and normalize bulk RNA-seq data
#' 
#' Subsets the genes to those that have more than a threshold number of counts per million (CPM) in a threshold number of samples. 
#' Then, normalizes the counts using the DESeq2 package.
#' @param df Data frame with raw counts. The rows should be genes, and columns should be samples. 
#' @param cell_thresh Minimum number of samples that need to have CPM above threshold (default=2)
#' @param count_thresh CPM threshold (default=0.5)
#' 
#' @return data frame with normalized counts
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

#' Filter and normalize scRNA-seq data
#' 
#' Filters the cells by the nUMI, nGene, log10GenesPerUMI, and mitoRatio threshold values,
#' gets the genes with non zero counts in X or more cells, where X=total cells * threshold proportion,
#' uses SCTransform to normalize the raw counts. Uses vst.flavor="v2",
#' and subsets the normalized counts table to the genes that pass threshold
#' 
#' @param df Data frame with raw counts. The rows should be genes, and columns should be samples. 
#' @param prop Threshold proportion of cells with non zero counts
#' @param nUMI_filt nUMI (total raw counts (nUMI) threshold
#' @param nGene_filt nGene (total number of non zero genes) threshold
#' @param log10_filt log10GenesPerUMI (non zero genes per total raw counts) threshold
#' @param mito_filt mitoRatio (mitochondrial gene counts per total raw counts) threshold
#' @param variable.features SCTransform subsets the dataset to the top n variable genes. This argument specifies n
#' @return data frame with normalized counts
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

#' Filter doublets from scRNA-seq data
#' 
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
  groups = ann_df %>% group_by(across(all_of(split.by))) %>% group_split() %>%
    lapply(FUN=function(x) {x %>% pull(!!sym(byvar))})
  df_list = lapply(groups, FUN=function(x) {
    df %>% select(one_of(x))
  })
  singlets = filter_doublets(df_list)
  print(head(singlets))
  df = df %>% select(one_of(singlets))
  return(df)
}
