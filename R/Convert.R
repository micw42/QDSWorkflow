#' @export
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

#' @export
convertSpecies = function(df, hm, old_id, new_id) {
  if (all(row.names(df)==seq(1, nrow(df)))) {
    gene_col = colnames(df)[1]
    df = df %>% distinct(!!sym(gene_col), .keep_all=T) %>%
      column_to_rownames(var=gene_col)
  }
  cols = colnames(hm)
  conf = cols[which(grepl("confidence", cols))]
  hom_type = cols[which(grepl("type", cols))]
  df = convert_species(df, hm=hm, hom_type=hom_type, 
                       conf=conf, old_id=old_id, new_id=new_id)
  return(df)
}