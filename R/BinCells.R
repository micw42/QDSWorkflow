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