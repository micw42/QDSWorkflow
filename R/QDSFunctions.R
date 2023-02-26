select = dplyr::select
rename = dplyr::rename
reduce = purrr::reduce
slice = dplyr::slice

normal_deseq = function(full_df, cond_df) {
  dds <- DESeqDataSetFromMatrix(countData = full_df, colData = cond_df, design = ~ 1)
  vsd <- DESeq2::varianceStabilizingTransformation(dds)
  vsd = as.data.frame(assay(vsd))
  return(vsd)
}

#V2 norm
v2_sct = function(df, variable.features=3000) {
  df = df %>%
    SCTransform(vst.flavor = "v2", variable.features.n=variable.features) 
  df_norm = df %>%
    GetAssayData(slot="scale.data")
  print("Got assay data.")
  rm(df)
  df_norm = df_norm %>% as.data.frame() 
  print("Converted to data frame.")
  return(df_norm)
}

seurat_norm = function(df) {
  df = df %>%
    NormalizeData()
  df = df %>%
    GetAssayData(slot="data") %>%
    as.data.frame()
  return(df)
}


filter_cells = function(df, nUMI_filt=500, nGene_filt=250, log10_filt=0.75, mito_filt=0.2) {
  filtered_seurat <- subset(x = df, 
                            subset= (nUMI >= nUMI_filt) & 
                              (nGene >= nGene_filt) & 
                              (log10GenesPerUMI > log10_filt) & 
                              (mitoRatio < mito_filt))
  return(filtered_seurat)
}

sc_preproc = function(df, prop=0.3, 
                      nUMI_filt=500, nGene_filt=250, log10_filt=0.75, 
                      mito_filt=0.2, variable.features=3000) {
  print("nGenes Before filter:")
  print(nrow(df))
  df = df %>% get_seurat_obj() %>% filter_cells(nUMI_filt=nUMI_filt, nGene_filt=nGene_filt,
                                                log10_filt=log10_filt, mito_filt=mito_filt)
  good_genes = df %>% get_good_genes(prop=prop)
  print(paste("number of good genes at prop =", prop, ":", length(good_genes)))
  df = df %>% v2_sct(variable.features=variable.features)
  print("nGenes after v2 sct:")
  print(nrow(df))
  df = df %>%
    rownames_to_column(var="Geneid") %>%
    filter(Geneid %in% good_genes)
  print("nGenes after good_genes filter:")
  print(nrow(df))
  return(df) }

get_seurat_obj = function(df, mt.pattern="^MT-") {
  df = df %>% 
    as.matrix() %>% CreateSeuratObject()
  df$log10GenesPerUMI <- log10(df$nFeature_RNA) / log10(df$nCount_RNA)
  df$mitoRatio <- PercentageFeatureSet(object = df, pattern = mt.pattern)
  df$mitoRatio <- df@meta.data$mitoRatio / 100
  metadata <- df@meta.data
  metadata <- metadata %>%
    dplyr::rename(sample = orig.ident,
                  nUMI = nCount_RNA,
                  nGene = nFeature_RNA)
  df@meta.data <- metadata
  return(df)
}

make_hist = function(obj, metric, grouping=NULL, ann_df=NULL) {
  meta = obj@meta.data
  if (!is.null(ann_df)) {
    meta = meta %>% rownames_to_column(var="sample_id") %>%
      inner_join(ann_df, by=c("sample_id"=colnames(ann_df)[1])) %>%
      column_to_rownames(var="sample_id")
  }
  p = meta %>%
    ggplot(aes_string(x=meta[[metric]], color = grouping, fill=grouping)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    scale_x_log10() +	
    xlab(metric) +
    ggtitle(paste(metric, "\ngrouped by", grouping))
  return(p)
}

extract_nzc <- function(fit, lm, family){
  coefs <- coef(fit, s=lm)
  if (family == "multinomial"){nzc <- unique(unlist(lapply(coefs, nzcs)))}
  else(nzc <- nzcs(coefs))
  
  return(nzc)
}

bulk_preproc = function(df, cell_thresh=2, count_thresh=0.5) {
  good_genes = df %>% cpm_filter(cell_thresh=cell_thresh, count_thresh=count_thresh) %>% rownames()
  cond_df = data.frame(row.names=colnames(df), Condition=rep(1, ncol(df)))
  df = (df %>% normal_deseq(cond_df = cond_df)) %>%
    rownames_to_column(var="Geneid") %>%
    filter(Geneid %in% good_genes)
  return(df)
}

get_good_genes = function(obj, prop=0.1) {
  counts = obj %>% GetAssayData()  
  thresh = (ncol(counts)*prop) %>% round()  
  nonzero <- counts > 0
  keep_genes <- Matrix::rowSums(nonzero) >= thresh
  good_genes = rownames(counts)[keep_genes]
  return(good_genes)
}


convert_species = function(full_df, hm, hom_type, conf, old_id, new_id) {
  hm = hm %>% filter(!!sym(hom_type)=="ortholog_one2one" & !!sym(conf)==1) %>%
    select(c(old_id, new_id)[1:2]) %>%
    distinct(!!sym(old_id), .keep_all=T) 		
  conv_df = full_df %>% rownames_to_column(var="Geneid") %>%
    inner_join(hm, by=c("Geneid"=old_id)) %>%
    select(-Geneid) %>%
    distinct(!!sym(new_id), .keep_all = T) %>%
    dplyr::rename(Geneid=!!sym(new_id)) %>%
    relocate(Geneid)
  return(conv_df)
}


# convert gene symbol to ensembl ID
convert_symbol = function(gene_map, full_df, conv_col="converted_alias", init_col="initial_alias") {
  gene_map = gene_map %>%
    select(one_of(c(conv_col, init_col)))
  gene_map = gene_map %>%
    rename("converted_alias"=conv_col, "initial_alias"=init_col) %>%	
    filter(converted_alias!="None") %>%
    distinct(initial_alias, .keep_all=T)
  
  full_df = full_df %>% rownames_to_column(var="Geneid") %>%
    inner_join(gene_map, by=c("Geneid"="initial_alias")) %>%
    select(-Geneid) %>%
    dplyr::rename(Geneid=converted_alias) %>%
    relocate(Geneid) %>%
    distinct(Geneid, .keep_all=T)
  
  return(full_df)
}

#df_list: named list of data frames to combine
run_ComBat = function(df_list, ann, batch_col, join_by="Geneid", cov_col=NULL, mean.only=F, ref.batch=NULL, int_counts=F) {
  col_names = df_list %>% lapply(FUN = function(x) {
    colnames(x)[-1]
  })
  full_df = df_list %>% reduce(inner_join, by = join_by) %>% column_to_rownames(var="Geneid")
  samples = colnames(full_df)
  ann = ann %>% slice(match(samples, sample_id))
  batch = ann %>% pull(batch_col)
  
  if (!is.null(cov_col)) {
    formula = paste0("~as.factor(", cov_col, ")")
    mod = model.matrix(as.formula(formula), data = ann)
  } else {
    mod=NULL
  }
  
  if (!int_counts) {
    combat_edata = ComBat(full_df, batch = batch, mean.only = mean.only, mod = mod, ref.batch=ref.batch)
  } else {
    full_df = as.matrix(full_df)
    combat_edata = ComBat_seq(full_df, batch = batch, covar_mod = mod)
  }
  out = lapply(col_names, FUN = function(x) {
    df = combat_edata %>% as.data.frame() %>% dplyr::select(one_of(x)) %>% t() %>% as.data.frame()
    return(df)
  })
  return(out)
}

nzcs <- function(coefs){
  nzc <- coefs@Dimnames[[1]][coefs@i + 1]
  nzc <- nzc[2:length(nzc)]
  
  return(nzc)
}


cpm_filter = function(df, cell_thresh=2, count_thresh=0.5) {
  cpm_df = cpm_norm(df)
  idx = rowSums( cpm_df >= count_thresh ) >= cell_thresh
  df = df[idx,]
  return(df)
}

cpm_norm = function(df) {
  cpm_df = df %>% scale(center = FALSE,
                        scale = colSums(df)/1e6)
  return(cpm_df)
}

y_weights <- function(xy){
  cts <- xy %>% dplyr::count(y)
  tot <- sum(cts$n)
  cts$Weight <- 1 - cts$n / tot; cts$n <- NULL; names(cts) <- c("y", "Weight")
  
  wdf <- inner_join(data.frame(y=xy$y), cts, by="y")
  weights <- wdf$Weight
  
  return(weights)
}

build_model <- function(xy, alpha, nfolds = 5, family="gaussian", type.measure="mse", weighted=T, standardize=T, keep=T){
  xy <- xy %>% arrange(y)
  x <- as.matrix(xy %>% select(-y))
  y <- as.matrix(xy %>% select(y))
  
  if (weighted){weights <- y_weights(xy)}
  else{weights <- NULL}
  
  fit <- cv.glmnet(x, y, nfolds=nfolds, family=family, type.measure=type.measure, weights=weights, alpha=alpha, standardize=standardize, keep=keep)
  
  return(fit)
}

alpha_test = function(x, y, n_alphas, nfolds=NULL, family="gaussian", type.measure="mse", 
                      err_range=1.05) {
  if (is.null(nfolds)) {
    nfolds=nrow(x)
  }
  x["y"] = y
  alphas = seq(0, 1, length.out=n_alphas)
  fit_list = list()
  err_vec = c()
  for (i in 1:length(alphas)) {
    alpha = alphas[i]
    fit = build_model(x, alpha=alpha, nfolds=nfolds, family=family, type.measure=type.measure)
    fit_list[[i]] = fit
    mce = min(fit$cvm)
    err_vec = c(err_vec, mce)
  }
  err_df = data.frame(Alpha=alphas, Errors=err_vec)
  min_err = min(err_df$Errors)
  alpha_range = err_df %>% filter(Errors<err_range*min_err) %>%
    pull(Alpha)
  opt_alpha = min(alpha_range)
  opt_idx = which(err_df$Alpha==opt_alpha)
  opt_model = fit_list[[opt_idx]]
  out = setNames(list(err_df, fit_list, opt_model, opt_alpha), 
                 c("err", "fit", "opt_model", "opt_alpha"))
  return(out)
}

best_preval = function(fit) {
  best <- fit$index["min",]
  pred = fit$fit.preval[, best]
  pred = exp(pred)/(1+exp(pred))
  return(pred)
}

cv_confusion = function(fit, ann_df, family="binomial", thresh=0.5) {
  best <- fit$index["min",]
  if (family=="binomial") {
    pred = fit$fit.preval[, best]
    pred = exp(pred)/(1+exp(pred))
    pred = prob_to_binom(pred, thresh=thresh)
    y = ann_df %>% slice(match(names(pred), sample_id)) %>%
      pull(Condition)	
    cnf = confusion.glmnet(pred, newy = y, family = family)
  } else if (family=="multinomial") {
    prevals = fit$fit.preval
    y = ann_df %>% slice(match(rownames(prevals), sample_id)) %>%
      pull(Condition)	
    cnf = confusion.glmnet(prevals, newy = y, family = "multinomial")
    cnf = cnf[[best]]
  }
  return(cnf)
}

cv_rates = function(fit, ann_df, pos="HR", family="binomial", thresh=0.5) {
  cnf_best <- cv_confusion(fit, ann_df, family=family, thresh=thresh)
  if (family=="multinomial") {
    return(cnf_best)
  } else if (family=="binomial") {		
    df = as.data.frame(cnf_best)
    cv_fnr_fpr = mat_to_err(df, pos=pos)
    out = setNames(list(cv_fnr_fpr, cnf_best), c("cv_rates", "cv_confusion"))
    return(out)
  }	
}

mat_to_err = function (err_df, pos) 
{
  tp = err_df %>% filter(Predicted == 1 & True == pos) %>% 
    pull(Freq)
  if (length(tp) == 0) {
    tp = 0
  }
  fp = err_df %>% filter(Predicted == 1 & True != pos) %>% 
    pull(Freq)
  if (length(fp) == 0) {
    fp = 0
  }
  tn = err_df %>% filter(Predicted != 1 & True != pos) %>% 
    pull(Freq)
  if (length(tn) == 0) {
    tn = 0
  }
  fn = err_df %>% filter(Predicted != 1 & True == pos) %>% 
    pull(Freq)
  if (length(fn) == 0) {
    fn = 0
  }
  fpr = fp/(fp + tn)
  fnr = fn/(fn + tp)
  mce = (fp + fn)/(fp + fn + tp + tn)
  out = setNames(list(mce, fpr, fnr), c("MCE", "FPR", "FNR"))
  return(out)
}

#' @export
make_grouped_hist = function(pred_df, ann_df, grouping1, grouping2, 
                             byvar=c("sample_id"="sample_id"), title=NULL, levels=NULL) {
  pred_df = inner_join(pred_df, ann_df, by=byvar) 
  if (!is.null(levels)) {
    pred_df[[grouping1]] = factor(pred_df[[grouping1]], ordered=T, levels=levels)
  }
  modes = pred_df %>% group_by(!!sym(grouping1), !!sym(grouping2)) %>%
    group_split() %>%
    lapply(FUN=function(x) {
      cond = x %>% pull(!!sym(grouping1)) %>% unique()
      mode_idx = which.max(density(x$QDS)$y)
      mode = density(x$QDS)$x[mode_idx]
      df = data.frame(Condition=cond, QDS_mode=mode)
      colnames(df)[1] = grouping1
      return(df)
    })
  modes_df = bind_rows(modes) %>% group_by(!!sym(grouping1)) %>% summarise(QDS_mod_med=median(QDS_mode))
  
  p = ggplot(pred_df, aes(x=QDS, fill=NULL, group=!!sym(grouping2)))+
    geom_density(adjust=1.5, alpha=.4) +
    geom_vline(data=modes_df, mapping=aes(xintercept=QDS_mod_med, color="red")) +
    facet_wrap(as.formula(paste("~", grouping1)), ncol=1)
  return(p)
}

convert_to_numeric = function(vec, pos=NULL) {
  num_vec = as.numeric(vec==pos)
  return(num_vec)
}

prob_to_binom = function(prob_vec, thresh=0.5) {
  pos_idx = which(prob_vec > thresh)
  neg_idx = which(prob_vec < thresh)
  prob_vec[pos_idx] = 1
  prob_vec[neg_idx] = 0
  return(prob_vec)
}

multinom_mce = function(pred, actual, lev=c("G", "Q", "S"), ignore=NULL) {
  if (!is.null(ignore)) {
    ignore_idx=(actual %in% ignore)
    print(ignore_idx)
    pred = pred[!ignore_idx]
    actual = actual[!ignore_idx]
  }	
  errs = list()
  for (i in 1:length(lev)) {
    class = lev[i]
    class_pred = pred %>% convert_to_numeric(pos=class)
    class_actual = actual %>% convert_to_numeric(pos=class)
    class_err = mce_fnr_fpr(actual = class_actual, pred = class_pred)
    errs[[i]] = class_err
  }
  n_incorrect = length(which(pred!=actual))
  overall_mce = n_incorrect/length(pred)
  names(errs) = lev
  mce_vec = lapply(errs, FUN=function(x) x$MCE) %>% unlist()
  fnr = lapply(errs, FUN=function(x) x$FNR) %>% unlist()
  fpr = lapply(errs, FUN=function(x) x$FPR) %>% unlist()
  err_df = data.frame(MCE=mce_vec, FNR=fnr, FPR=fpr)
  out = setNames(list(err_df, overall_mce), c("err_df", "overall_mce")) 
  return(out) 
}



wrap_box = function(df, ann_df, byvar=c("sample_id"="sample_id"), x, y, wrap, color=NULL, title=NULL) {
  df = inner_join(df, ann_df, by=byvar)
  p = ggplot(df, aes(x = !!sym(x), y = !!sym(y))) + 
    geom_boxplot(outlier.shape = NA)
  if (!is.null(color)) {
    p = p + geom_jitter(aes(col=!!sym(color)), shape = 16, height = 0, width = 0.2)
  } else {
    p = p + geom_jitter(shape = 16, height = 0, width = 0.2)	 
  }
  p = p + facet_wrap(as.formula(paste("~", wrap))) +
    ggtitle(title)
  return(p)
}

#' @export
make_boxplot = function(df, ann_df, x, y="QDS", byvar=c("sample_id"="sample_id"), color=NULL, title=NULL, 
                        levels=NULL) {
  df = inner_join(df, ann_df, by=byvar) 
  if (!is.null(levels)) {
    df[[x]] = factor(df[[x]], ordered=T, levels=levels)
  }
  p = ggplot(df, aes_string(x = x, y = y)) + 
    geom_boxplot(outlier.shape = NA) + geom_jitter(shape = 16, 
                                                   height = 0, width = 0.2, aes_string(col=color)) + ggtitle(title)
  return(p)
}

max_fnr_fpr = function(actual, pred) {
  pos = unique(actual)[1] %>% as.character()
  actual = actual %>% convert_to_numeric(pos=pos)
  pred = pred %>% convert_to_numeric(pos=pos)
  err = mce_fnr_fpr(pred=pred, actual=actual)
  m = max(err$FNR, err$FPR)
  return(m)
}

summ_func = function(data, lev, model) {
  actual = data$obs
  pred = data$pred
  metric = max_fnr_fpr(actual=actual, pred=pred)
  return(metric)
}

w_train_test = function(ann_df, grouping=c("Dataset", "Condition"), k_folds=5) {
  ann_list = ann_df %>% group_by_at(grouping) %>%
    group_split()
  fold_list = list()
  for (i in 1:length(ann_list)) {
    df = ann_list[[i]]
    nrow_i       <- nrow(df)
    n_per_fold_i <- ceiling(nrow_i / k_folds)
    df$fold = sample(rep(1:k_folds, n_per_fold_i), nrow_i, replace = FALSE)
    fold_list[[i]] = df
  }
  folds = fold_list %>% bind_rows() %>%
    group_by(fold) %>% group_split() %>%
    lapply(FUN=function(x) x$sample_id)
  return(folds)
}

bin_cells = function (df, n_cells = 60, bin_function="sum") 
{
  n_bin = round(nrow(df)/n_cells)
  df$bin = sample(seq(1, n_bin), nrow(df), replace = T)
  binned_df = df %>% group_by(bin) %>% summarise_if(is.numeric, 
                                                    bin_function) %>% select(-bin)
  print("Done binning.")
  return(binned_df)
}

format_bin = function(x, group_col, n_cells, bin_function) {
  sample_origin = (x[[group_col]])[1]
  x = x %>% bin_cells(n_cells = n_cells, bin_function=bin_function)
  ids = paste(sample_origin, seq(1, nrow(x)), sep = "_")
  x = x %>% mutate(sample_id = ids)
  x[group_col] = sample_origin
  return(x)
}

bin_wrapper = function (df, ann_df, join_col, group_col, n_cells = 60, bin_function="sum",
                        parallel=F, par_idx=NULL) 
{
  df_list = inner_join(df, ann_df, by = join_col) %>% group_by(!!sym(group_col)) %>% 
    group_split()
  if (!parallel) {       
    binned_df_list = lapply(df_list, FUN = function(x) {
      x = format_bin(x, group_col, n_cells, bin_function)
      return(x)
    }) 
    binned_df = bind_rows(binned_df_list)
  } else {
    x = df_list[[par_idx]]
    binned_df = format_bin(x, group_col, n_cells, bin_function)
  }    
  counts = binned_df %>% select(-c(group_col)) %>% column_to_rownames(var = "sample_id") %>% 
    t() %>% as.data.frame() %>% rownames_to_column(var = "Geneid")
  ann_df = binned_df %>% select(c("sample_id", group_col))
  out = setNames(list(counts, ann_df), c("counts", "ann_df"))
  return(out)
}

summarise_group = function(df, ann_df, group_col, byvar="Index", bin_function="median") {
  df = inner_join(df, ann_df, by = byvar) %>%
    group_by(!!sym(group_col)) %>%
    summarise_if(is.numeric, bin_function)
  return(df)
}

filter_doublets = function(obj, ann_df, split.by="Sample", 
                           byvar="Index") {
  meta = obj@meta.data %>% rownames_to_column(var="Index") %>%
    inner_join(ann_df, by=byvar) %>%
    column_to_rownames(var="Index")
  print(paste("N total:", nrow(meta)))
  obj@meta.data = meta
  obj.split <- SplitObject(obj, split.by = split.by) 
  singlet_list = list()
  for (i in 1:length(obj.split)) {
    df = obj.split[[i]]        
    if (ncol(df)<50) {
      meta_df = df@meta.data
      singlets <- meta_df %>% rownames()
      singlet_list[[i]] = singlets
    } else {
      df = df %>% NormalizeData() %>%
        FindVariableFeatures() %>%
        ScaleData() %>%
        RunPCA(nfeatures.print = 10, npcs=min(50, (ncol(df)-1)))
      
      # Find significant PCs
      stdv <- df[["pca"]]@stdev
      sum.stdv <- sum(df[["pca"]]@stdev)
      percent.stdv <- (stdv / sum.stdv) * 100
      cumulative <- cumsum(percent.stdv)
      co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
      co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                           percent.stdv[2:length(percent.stdv)]) > 0.1), 
                  decreasing = T)[1] + 1
      min.pc <- min(co1, co2)
      
      df <- df %>% RunUMAP(dims = 1:min.pc, n_neighbors=min(30, (ncol(df)-1))) %>%
        FindNeighbors(dims = 1:min.pc) %>%           
        FindClusters(resolution = 0.1)
      sweep.list <- paramSweep_v3(df, PCs = 1:min.pc, num.cores = detectCores() - 1)
      sweep.stats <- summarizeSweep(sweep.list)
      bcmvn <- find.pK(sweep.stats)
      print("Done with sweep.")
      
      bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
      optimal.pk <- bcmvn.max$pK
      optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
      print("Found optimal pk.")
      
      annotations <- df@meta.data$seurat_clusters
      homotypic.prop <- modelHomotypic(annotations) 
      nExp.poi <- round(optimal.pk * nrow(df@meta.data)) 
      nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
      print("Found nExp.poi.")
      
      # run DoubletFinder
      df <- doubletFinder_v3(seu = df, 
                             PCs = 1:min.pc, 
                             pK = optimal.pk,
                             nExp = nExp.poi.adj)
      
      meta_df = df@meta.data
      colnames(meta_df)[grepl("DF.classification", colnames(meta_df))] = "doublet_finder"
      singlets <- meta_df %>% filter(doublet_finder == "Singlet") %>% rownames()
      singlet_list[[i]] = singlets
    }
    
  }
  singlets <- unlist(singlet_list)
  print(paste("N singlets:", length(singlets)))
  return(singlets)
}


make_time_series = function(df, ann_df, grouping, sum_var,
                            byvar = "sample_id",
                            summ_func="mean", title=NULL, group=1) {
  df = inner_join(df, ann_df, by = byvar) %>%
    group_by(!!sym(grouping)) %>%
    summarise_if(is.numeric, summ_func)
  print(df)
  p = ggplot(data=df, aes_string(x=grouping, y=sum_var, group=group)) +
    geom_line()+
    geom_point()+
    ggtitle(title)
  return(p)
}

format_cols = function(df) {
  if (all(row.names(df)==seq(1, nrow(df)))) {
    gene_col = colnames(df)[1]
    df = df %>% distinct(!!sym(gene_col), .keep_all=T) %>%
      column_to_rownames(var=gene_col)
  }
  return(df)
}