#' Gene symbol conversion
#' 
#' Converts the gene symbols in the counts table to Ensembl IDs, or Ensembl IDs into gene symbols
#' 
#' @param df Data frame with raw counts. The rows should be genes, and columns should be samples. 
#' @param gene_map Gene table from gProfiler, loaded as data frame
#' @param conv_col Name of column in gene_map with converted symbols/IDs
#' @param init_col Name of column in gene_map with original symbols/IDs
#' 
#' @return New data frame with converted genes in the column called "Geneid"
#' @export
convertSymbol = function(df, gene_map, conv_col="converted_alias", init_col="initial_alias") {
  if (all(row.names(df)==seq(1, nrow(df)))) {
    gene_col = colnames(df)[1]
    df = df %>% distinct(!!sym(gene_col), .keep_all=T) %>%
      column_to_rownames(var=gene_col)
  }
  df = convert_symbol(gene_map, df, conv_col=conv_col, init_col=init_col)
  return(df)
}

#' Convert gene IDs between species
#' 
#' Convert the IDs in the counts table to the IDs of the new species (one-to-one, high confidence homologs only)
#' 
#' @param df Data frame with raw counts. The rows should be genes, and columns should be samples
#' @param hm Homology table from Ensembl website, loaded as data frame
#' @param old_id Name of column in homology table containing the original species IDs
#' @param new_id Name of column in homology table containing the new species IDs
#' @return New data frame with converted genes in the column called "Geneid"
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