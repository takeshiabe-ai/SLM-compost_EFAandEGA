pkgs <- c("lavaan","readxl","stats")
invisible(lapply(pkgs, function(p){
  if(!requireNamespace(p, quietly=TRUE)) stop("Missing package: ", p)
}))
library(lavaan); library(readxl)

dat <- read.csv(data_file)
dat_num <- numify_df(dat)

syn_tbl <- readxl::read_excel(syntax_xlsx, sheet=sheet_name)
syn_tbl <- syn_tbl[match(need_models, syn_tbl$model), ]

extract_factor_item_counts <- function(model_syntax){
  syn <- clean_syntax(model_syntax)
  lines <- unlist(strsplit(syn, "\n", fixed=TRUE))
  lines <- trimws(lines)
  lines <- lines[lines != ""]
  lines <- lines[!grepl("^#", lines)]
  meas  <- lines[grepl("=~", lines, fixed=TRUE)]
  if(length(meas) == 0){
    return(data.frame(factor=character(0), n_items=integer(0)))
  }
  out <- lapply(meas, function(line){
    lhs <- trimws(sub("=~.*", "", line))
    rhs <- trimws(sub("^.*=~", "", line))
    parts <- trimws(unlist(strsplit(rhs, "\\+")))
    parts <- parts[parts != ""]
    items <- trimws(sub(".*\\*", "", parts))
    items <- items[items != ""]
    data.frame(factor=lhs, n_items=length(unique(items)), stringsAsFactors=FALSE)
  })
  df <- do.call(rbind, out)
  aggregate(n_items ~ factor, df, sum)
}

items_per_factor_summary <- function(model_syntax){
  fc <- extract_factor_item_counts(model_syntax)
  if(nrow(fc) == 0) return(list(items_mean=NA_real_, items_se=NA_real_))
  n <- fc$n_items
  list(
    items_mean = mean(n),
    items_se   = if(length(n) > 1) sd(n) / sqrt(length(n)) else 0
  )
}

latent_cor_summary <- function(fit){
  lv <- lavaan::lavNames(fit, "lv")
  if(length(lv) < 2) return(list(cor_mean=NA_real_, cor_sd=NA_real_))

  pe <- lavaan::parameterEstimates(fit, standardized=FALSE)
  cov <- pe[pe$op=="~~" & pe$lhs!=pe$rhs & pe$lhs %in% lv & pe$rhs %in% lv, , drop=FALSE]
  if(nrow(cov) == 0) return(list(cor_mean=NA_real_, cor_sd=NA_real_))

  cov$key <- apply(cbind(cov$lhs, cov$rhs), 1, function(x) paste(sort(x), collapse="~~"))
  cov <- cov[!duplicated(cov$key), , drop=FALSE]

  list(
    cor_mean = mean(cov$est, na.rm=TRUE),
    cor_sd   = stats::sd(cov$est, na.rm=TRUE)
  )
}

run_cfa_metrics_one <- function(model_name, model_syntax, dat_num, min_n_cat=5){
  syn <- clean_syntax(model_syntax)
  items <- extract_items_from_syntax(syn)

  it_sum <- items_per_factor_summary(syn)

  out <- data.frame(
    model = model_name,
    N = NA_integer_,
    items_per_factor_mean = it_sum$items_mean,
    items_per_factor_se   = it_sum$items_se,
    lv_cor_mean = NA_real_,
    lv_cor_sd   = NA_real_,
    parameterization = parameterization_used,
    note = NA_character_,
    stringsAsFactors = FALSE
  )

  miss <- setdiff(items, colnames(dat_num))
  if(length(miss) > 0){
    out$note <- paste0("ITEM_NOT_FOUND: ", paste(miss, collapse=", "))
    return(out)
  }

  d <- dat_num[, items, drop=FALSE]
  cc <- complete.cases(d)
  d <- d[cc, , drop=FALSE]
  out$N <- nrow(d)
  if(nrow(d) < 10){
    out$note <- "TOO_FEW_COMPLETE_CASES"
    return(out)
  }

  d2 <- as.data.frame(lapply(d, collapse_sparse_ordinal, min_n=min_n_cat))
  d_ord <- as.data.frame(lapply(d2, make_ordered_observed_levels))

  bad <- names(which(sapply(d_ord, function(x) length(unique(x[!is.na(x)])) < 2)))
  if(length(bad) > 0){
    out$note <- paste0("SINGLE_CATEGORY_ITEM: ", paste(bad, collapse=", "))
    return(out)
  }

  fit <- tryCatch(
    lavaan::cfa(
      model = syn,
      data  = d_ord,
      ordered = colnames(d_ord),
      estimator = "WLSMV",
      parameterization = parameterization_used,
      std.lv = TRUE,
      missing = "listwise",
      control = list(iter.max=50000)
    ),
    error=function(e) e
  )

  if(inherits(fit, "error")){
    out$note <- paste0("CFA_ERROR: ", fit$message)
    return(out)
  }
  if(!isTRUE(lavaan::inspect(fit, "converged"))){
    out$note <- "NOT_CONVERGED"
    return(out)
  }

  cor_sum <- latent_cor_summary(fit)
  out$lv_cor_mean <- cor_sum$cor_mean
  out$lv_cor_sd   <- cor_sum$cor_sd
  out$note <- "OK"
  out
}

res_list <- lapply(seq_len(nrow(syn_tbl)), function(i){
  run_cfa_metrics_one(syn_tbl$model[i], syn_tbl$syntax[i], dat_num, min_n_cat=5)
})

res_tbl <- do.call(rbind, res_list)

out_csv <- file.path(out_sum, paste0("CFA_metrics_6models_WLSMV_", parameterization_used, ".csv"))
write.csv(res_tbl, out_csv, row.names=FALSE)

print(res_tbl)
cat("05_cfa_metrics_6models: saved metrics table.\n")
