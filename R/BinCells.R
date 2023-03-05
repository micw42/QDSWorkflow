#' Bin scRNA-seq counts
#' 
#' Groups the cells, partitions each group into random bins, and calculates a summary statistic for the counts in each bin.
#' @param df Data frame of normalized, filtered counts. The rows should be genes, and columns should be samples. 
#' @param ann_df Data frame containing the group of each cell. One column should contain cell IDs, rest of the columns should contain corresponding groups of the cell
#' @param group_col Name of column in ann_df containing the group of each cell
#' @param join_col Name of column in ann_df containing the cell IDs
#' @param n_bin Number of random bins to partition in each group
#' @param bin_function The function to use when calculating the summary statistic
#' @return Data frame with 1st column called "Geneid" containing gene names, rest of the columns containing the summary statistic for each bin. The bins are named by a random cell ID in the bin
#' @export
BinCells = function (df, ann_df, group_col, join_col="sample_id", n_bin = 10, bin_function="median") 
{
  df = df %>% format_cols() %>% t() %>% as.data.frame() %>% rownames_to_column(var="sample_id")
  df_list = inner_join(df, ann_df, by = c("sample_id"=join_col)) %>% group_by(!!sym(group_col)) %>% 
    group_split()
  binned_df_list = lapply(df_list, FUN = function(x) {
    x = x %>% bin_cells(n_bin=n_bin, bin_function=bin_function)
    return(x)
  }) 
  binned_df = bind_rows(binned_df_list)
  counts = binned_df %>% column_to_rownames(var = "sample_id") %>% 
    t() %>% as.data.frame() %>% rownames_to_column(var = "Geneid")
  return(counts)
}

#' Bootstrap scRNA-seq cells
#' 
#' Samples n cells with replacement from each group and taking the sum of their normalized counts
#' @param df Data frame of normalized, filtered counts. The rows should be genes, and columns should be samples. 
#' @param ann_df Data frame containing the group of each cell. One column should contain cell IDs, rest of the columns should contain corresponding groups of the cell
#' @param group_col Name of column(s) in ann_df to use when grouping the cells. Can be a string (for one column) or vector (for multiple columns).
#' @param join_col Name of column in ann_df containing the cell IDs
#' @param n Number of cells to sample with replacement in each bootstrap sample
#' @param m Number of bootstrap samples to generate in each group
#' @return A named list
#' \itemize{
#'   \item samples - Data frame of bootstrapped samples, with first column as gene ID
#'   \item ann - New ann_df, corresponding to the new sample IDs. Use this when ann_df is needed downstream.
#' }
#' @export
BootstrapSamples = function(df, ann_df, group_col, join_col="sample_id", n=50, m=100) {
  df = df %>% format_cols() %>%
    t() %>% as.data.frame() %>% rownames_to_column(var="sample_id")
  ann_df = ann_df %>% dplyr::select(one_of(c(join_col, group_col)))
  df_list = inner_join(df, ann_df, by = c("sample_id"=join_col)) %>% group_by(across(all_of(group_col))) %>% 
    group_split()
  bootstrap_df_list = list()
  map_list = list()
  for (i in 1:length(df_list)) {
    sub_df = df_list[[i]]
    cell = (sub_df %>% pull(sample_id))[1]
    cells = rep(cell, m)
    subset = sub_df %>% select(-one_of(c(group_col))) %>%
      column_to_rownames(var="sample_id") %>%
      bootstrap_cells(n=n, m=m)
    new_colnames = paste("bootstrap", i, seq(1, m), sep="_")
    colnames(subset) = new_colnames
    map_rows = data.frame(sample_id=new_colnames, old_id=cells)
    bootstrap_df_list[[i]] = subset
    map_list[[i]] = map_rows
  }
  
  bootstrap_df = bind_cols(bootstrap_df_list)
  map_df = bind_rows(map_list)
  new_ann = inner_join(map_df, ann_df, by=c("old_id"=join_col))
  
  out = setNames(list(bootstrap_df, new_ann), c("samples", "ann"))
  return(out)
}