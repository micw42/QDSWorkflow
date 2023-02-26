#' @export
QDS = function(train_df, test_df, ann_df, y_col) {
  train_df = train_df %>% format_cols()
  test_df = test_df %>% format_cols()
  train_ids = colnames(train_df)
  test_ids = colnames(test_df)
  train_df = train_df %>% rownames_to_column(var="Geneid")
  test_df = test_df %>% rownames_to_column(var="Geneid")
  
  combat_ann = data.frame(sample_id=c(train_ids, test_ids), Dataset=rep(c("train", "test"),
                                                                        times=c(length(train_ids), 
                                                                                length(test_ids))))
  out = run_ComBat(list(train_df, test_df), ann = combat_ann, batch_col = "Dataset")
  
  train_df = out[[1]]
  test_df = out[[2]]
  
  id_col = colnames(ann_df)[1]
  y = ann_df %>% 
    slice(match(train_ids, !!sym(id_col))) %>%
    pull(!!sym(y_col))
  out = alpha_test(train_df, y, n_alphas=20)
  model = out$opt_model
  pred_df = data.frame(sample_id=rownames(test_df),
                       pred=predict(model, newx=as.matrix(test_df),
                                    s=model$lambda.min, type="response")) %>%
    rename("QDS"="s1")
  return(pred_df)
}

#' @export
PrebuiltQDSModel = function(test_df) {
  train_df = readRDS(system.file("extdata", "REF_final_norm.rds", package="QDSWorkflow")) %>% format_cols()
  test_df = format_cols(test_df)
  model = readRDS(system.file("extdata", "ref_model.rds", package="QDSWorkflow"))
  train_ids = colnames(train_df)
  test_ids = colnames(test_df)
  train_df = train_df %>% rownames_to_column(var="Geneid")
  test_df = test_df %>% rownames_to_column(var="Geneid") 
  test_df = left_join(train_df, test_df, by="Geneid") %>%
    replace(is.na(.), 0) %>%
    dplyr::select(-one_of(train_ids))
  
  feat = extract_nzc(model, lm=model$lambda.min, family="gaussian")
  sub_df = test_df %>% slice(match(Geneid, feat)) %>% column_to_rownames(var="Geneid")
  print(sub_df)
  n_zero = mean(colSums((sub_df==0))/nrow(sub_df)) 
  print(n_zero)
  if (n_zero > 0.5) {
    warning("Samples are missing >20% selected features on average. We recommend building a new QDS model.")
  }
  combat_ann = data.frame(sample_id=c(train_ids, test_ids), Dataset=rep(c("train", "test"),
                                                                        times=c(length(train_ids), 
                                                                                length(test_ids))))
  out = run_ComBat(list(train_df, test_df), ann = combat_ann, batch_col = "Dataset", ref.batch = "train")
  test_df = out[[2]]
  pred_df = data.frame(sample_id=rownames(test_df),
                       pred=predict(model, newx=as.matrix(test_df),
                                    s=model$lambda.min, type="response")) %>%
    rename("QDS"="s1")
  return(pred_df)
}