convertSymbol = function(df, gene_map, conv_col="converted_alias", init_col="initial_alias") {
  if (all(row.names(df)==seq(1, nrow(df)))) {
    gene_col = colnames(df)[1]
    df = df %>% distinct(!!sym(gene_col), .keep_all=T) %>%
      column_to_rownames(var=gene_col)
  }
  rownames(df) = toupper(rownames(df))
  df = convert_symbol(gene_map, df, conv_col=conv_col, init_col=init_col)
  return(df)
}

