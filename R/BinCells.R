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