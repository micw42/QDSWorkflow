#' Train a model and predict quiescence depth
#' 
#' Merges the train and test datasets and removes batch effects between them. 
#' Then, builds an elastic net regression model on the train data. In elastic net, the alpha parameter controls how many coefficients (genes) are in the model: 
#' lower alpha gives models with more coefficients, higher alpha gives models with fewer coefficients. 
#' Tests a range of alpha values and finds the cross-validation error of each one. T
#' To capture a broader signature of quiescence depth, we use the minimum alpha with cross-validation error within 5% of the minimum error. 
#' 
#' @param test_df  test data (normalized counts). The rows should be genes, and columns should be samples.
#' @param train_df  train data (normalized counts). The rows should be genes, and columns should be samples.
#' @param ann_df data frame with train data labels. First column contains train sample IDs and 2nd column contains corresponding y values. 
#' @param y_col Name of column in ann_df that contains y values
#' @param prebuilt Set to TRUE if using the prebuilt rat embryonic fibroblast model to predict on the test data.  If prebuilt=TRUE, test_df is the only required argument.
#' @return Data frame with predicted quiescence depth value for each test sample
#' @export
QDS = function(test_df, train_df=NULL, ann_df=NULL, y_col=NULL, prebuilt=F) {
  if (prebuilt) {
    pred_df = PrebuiltQDSModel(test_df)
    return(pred_df)
  }
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
  out = suppressWarnings(alpha_test(train_df, y, n_alphas=20))
  model = out$opt_model
  pred_df = data.frame(sample_id=rownames(test_df),
                       pred=predict(model, newx=as.matrix(test_df),
                                    s=model$lambda.min, type="response")) %>%
    rename("QDS"="s1")
  return(list(out, pred_df))
}


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
  n_zero = mean(colSums((sub_df==0))/nrow(sub_df)) 
  if (n_zero > 0.6) {
    warning("Samples are missing >60% selected features on average. We recommend building a new QDS model.")
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