pkgs <- c("lavaan","readxl")
invisible(lapply(pkgs, function(p){
  if(!requireNamespace(p, quietly=TRUE)) stop("Missing package: ", p)
}))
library(lavaan); library(readxl)

dat <- read.csv(data_file)
dat_num <- numify_df(dat)

syn_tbl <- readxl::read_excel(syntax_xlsx, sheet=sheet_name)

if(!("model" %in% names(syn_tbl)) || !("syntax" %in% names(syn_tbl))){
  stop("Excel must contain columns: 'model' and 'syntax'.")
}
if(!all(need_models %in% syn_tbl$model)){
  stop("Missing required models in Excel: ",
       paste(setdiff(need_models, syn_tbl$model), collapse=", "))
}
syn_tbl <- syn_tbl[match(need_models, syn_tbl$model), ]

run_cfa_one <- function(model_name, model_syntax, dat_num, min_n_cat=5){
  syn <- clean_syntax(model_syntax)
  items <- extract_items_from_syntax(syn)

  out_fit <- data.frame(
    model=model_name,
    N=NA_integer_,
    n_items=length(items),
    chisq=NA_real_, df=NA_real_, pvalue=NA_real_,
    CFI=NA_real_, TLI=NA_real_, RMSEA=NA_real_, SRMR=NA_real_,
    parameterization = parameterization_used,
    note=NA_character_,
    stringsAsFactors=FALSE
  )

  miss <- setdiff(items, colnames(dat_num))
  if(length(miss) > 0){
    out_fit$note <- paste0("ITEM_NOT_FOUND: ", paste(miss, collapse=", "))
    return(out_fit)
  }

  d <- dat_num[, items, drop=FALSE]
  cc <- complete.cases(d)
  d <- d[cc, , drop=FALSE]
  out_fit$N <- nrow(d)

  if(nrow(d) < 10){
    out_fit$note <- "TOO_FEW_COMPLETE_CASES"
    return(out_fit)
  }

  d2 <- as.data.frame(lapply(d, collapse_sparse_ordinal, min_n=min_n_cat))
  d_ord <- as.data.frame(lapply(d2, make_ordered_observed_levels))

  bad <- names(which(sapply(d_ord, function(x) length(unique(x[!is.na(x)])) < 2)))
  if(length(bad) > 0){
    out_fit$note <- paste0("SINGLE_CATEGORY_ITEM: ", paste(bad, collapse=", "))
    return(out_fit)
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
    out_fit$note <- paste0("CFA_ERROR: ", fit$message)
    return(out_fit)
  }
  if(!isTRUE(lavaan::inspect(fit, "converged"))){
    out_fit$note <- "NOT_CONVERGED"
    return(out_fit)
  }

  fm <- lavaan::fitMeasures(fit, c("chisq","df","pvalue","cfi","tli","rmsea","srmr"))
  out_fit$chisq  <- as.numeric(fm["chisq"])
  out_fit$df     <- as.numeric(fm["df"])
  out_fit$pvalue <- as.numeric(fm["pvalue"])
  out_fit$CFI    <- as.numeric(fm["cfi"])
  out_fit$TLI    <- as.numeric(fm["tli"])
  out_fit$RMSEA  <- as.numeric(fm["rmsea"])
  out_fit$SRMR   <- as.numeric(fm["srmr"])
  out_fit$note <- paste0("OK(WLSMV/", parameterization_used, ")")
  out_fit
}

fit_list <- lapply(seq_len(nrow(syn_tbl)), function(i){
  run_cfa_one(syn_tbl$model[i], syn_tbl$syntax[i], dat_num, min_n_cat=5)
})

fit_tbl <- do.call(rbind, fit_list)

out_csv <- file.path(out_cfa, paste0("CFA_fit_6models_WLSMV_", parameterization_used, ".csv"))
write.csv(fit_tbl, out_csv, row.names=FALSE)

print(fit_tbl)
cat("04_cfa_6models: saved fit table.\n")
