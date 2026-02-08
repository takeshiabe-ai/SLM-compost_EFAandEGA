pkgs <- c("EGAnet","lavaan","psych","qgraph","igraph","grDevices")
invisible(lapply(pkgs, function(p){
  if(!requireNamespace(p, quietly=TRUE)) stop("Missing package: ", p)
}))
library(EGAnet); library(lavaan); library(psych); library(qgraph); library(igraph)

# -----------------------------
# Helpers (EGA-specific)
# -----------------------------
to_ordered_no_empty <- function(x){
  x <- numify(x)
  lev <- sort(unique(x[!is.na(x)]))
  if(length(lev) < 2) return(NULL)  # drop constants
  ordered(x, levels=lev)
}

ensure_wc_named <- function(wc, items_use){
  wc <- as.integer(wc)
  if(is.null(names(wc)) || any(names(wc)=="")){
    if(length(wc) != length(items_use)){
      stop("wc length != items_use length.")
    }
    names(wc) <- items_use
  }
  wc
}

relabel_wc_1k <- function(wc_named){
  u <- sort(unique(as.integer(wc_named)))
  map <- setNames(seq_along(u), u)
  wc2 <- unname(map[as.character(as.integer(wc_named))])
  names(wc2) <- names(wc_named)
  wc2
}

make_comm_palette <- function(K){
  pal <- grDevices::palette.colors(K, palette = "Set3")
  names(pal) <- as.character(seq_len(K))
  pal
}

reassign_singletons <- function(wc_named, dat_num){
  cl <- split(names(wc_named), wc_named)
  single_ids <- names(cl)[sapply(cl, length) == 1]
  if(length(single_ids) == 0) return(list(wc=wc_named, note="no_singletons"))

  R <- suppressWarnings(cor(dat_num[, names(wc_named), drop=FALSE], use="pairwise.complete.obs"))
  for(cid in single_ids){
    item <- cl[[cid]][1]
    other_ids <- setdiff(names(cl), cid)
    scores <- sapply(other_ids, function(oid){
      members <- cl[[oid]]
      mean(abs(R[item, members]), na.rm=TRUE)
    })
    best <- other_ids[which.max(scores)]
    wc_named[item] <- as.integer(best)
    cl <- split(names(wc_named), wc_named)
  }
  list(wc=wc_named, note=paste0("reassigned_singletons=", length(single_ids)))
}

build_cfa_syntax <- function(wc_named, prefix="E"){
  sp <- split(names(wc_named), wc_named)
  lines <- character(0)
  k <- 0
  for(items in sp){
    k <- k + 1
    f <- paste0(prefix, k)
    if(length(items) == 1){
      next
    } else if(length(items) == 2){
      lab <- paste0("l", k)
      lines <- c(lines, paste0(f, " =~ ", lab, "*", items[1], " + ", lab, "*", items[2]))
    } else {
      lines <- c(lines, paste0(f, " =~ ", paste(items, collapse=" + ")))
    }
  }
  paste(lines, collapse="\n")
}

plot_rep_comm <- function(rep_comm, pal, main){
  if(is.null(rep_comm) || all(!is.finite(rep_comm))){
    plot.new(); text(0.5,0.5,"Replication not available", cex=1.1)
    mtext(main, side=3, line=1)
    return(invisible(NULL))
  }
  rep_comm <- rep_comm[order(as.integer(names(rep_comm)))]
  cols <- pal[names(rep_comm)]
  barplot(rep_comm, col=cols, ylim=c(0,1), las=1,
          ylab="Mean item replication (within community)",
          xlab="Community", main=main,
          names.arg=paste0("C", names(rep_comm)))
  abline(h=0.50, lty=2)
  legend("topright", legend=paste0("C", names(rep_comm)), fill=cols, bty="n", cex=0.9)
  invisible(NULL)
}

# -----------------------------
# Load & prep
# -----------------------------
dat <- read.csv(data_file)

items_num <- dat[, vars_full, drop=FALSE]
items_num <- numify_df(items_num)
items_num <- na.omit(items_num)
N_use <- nrow(items_num)

ord_list <- lapply(items_num, to_ordered_no_empty)
keep <- !sapply(ord_list, is.null)
dropped_const <- names(ord_list)[!keep]
dat_ord <- as.data.frame(ord_list[keep])
items_use <- names(dat_ord)

if(length(items_use) < 6) stop("Too few usable items after dropping constants.")

# -----------------------------
# 1) EGA
# -----------------------------
ega_res <- EGAnet::EGA(dat_ord, model="glasso", corr="auto", plot.EGA=FALSE)

wc0 <- ensure_wc_named(ega_res$wc, items_use)
wc_ega <- relabel_wc_1k(wc0)[items_use]
K <- length(unique(wc_ega))

pal <- make_comm_palette(K)
groups_list <- split(seq_along(items_use), wc_ega)
names(groups_list) <- paste0("C", seq_len(K))

# -----------------------------
# 2) bootEGA + itemStability
# -----------------------------
boot_res <- EGAnet::bootEGA(
  dat_ord,
  model="glasso", corr="auto",
  iter=boot_iter, type="resampling",
  seed=123,
  plot.itemStability=FALSE,
  typicalStructure=FALSE,
  plot.typicalStructure=FALSE,
  verbose=FALSE
)

IS <- EGAnet::itemStability(boot_res, IS.plot = FALSE)

rep_item <- NULL
if(!is.null(IS$item.stability) && !is.null(IS$item.stability$empirical.dimensions)){
  rep_item <- IS$item.stability$empirical.dimensions
  rep_item <- rep_item[items_use]
}

rep_comm <- NULL
if(!is.null(rep_item)){
  rep_comm <- tapply(rep_item, wc_ega, mean, na.rm=TRUE)
  rep_comm <- rep_comm[as.character(seq_len(K))]
  names(rep_comm) <- as.character(seq_len(K))
}

# Item replication PNG
png_file <- file.path(out_ega, "Replication_bootEGA_items.png")
png(png_file, width=1600, height=1200, res=200)
if(is.null(rep_item) || all(!is.finite(rep_item))){
  plot.new(); text(0.5,0.5,"Item replication not available", cex=1.1)
} else {
  rep_item2 <- sort(rep_item, decreasing = FALSE)
  cols <- pal[as.character(wc_ega[names(rep_item2)])]
  barplot(rep_item2, col=cols, las=2, ylim=c(0,1),
          ylab="Proportion replicated (item in its dimension)",
          main=paste0("Item Replication (bootEGA)  B=", boot_iter))
  abline(h=0.50, lty=2)
  legend("topleft", legend=paste0("Community ", 1:K), fill=pal[as.character(1:K)], bty="n", cex=0.9)
}
dev.off()

# -----------------------------
# 3) CFA (WLSMV / delta or theta)
# -----------------------------
wc_fix_raw <- reassign_singletons(wc0, items_num)$wc
wc_fix <- relabel_wc_1k(wc_fix_raw)
model_syn <- build_cfa_syntax(wc_fix, prefix="E")

fit_row <- data.frame(
  N = N_use,
  n_items = length(items_use),
  n_clusters_EGA = K,
  n_clusters_for_CFA = length(unique(wc_fix)),
  chisq=NA_real_, df=NA_real_, pvalue=NA_real_,
  CFI=NA_real_, TLI=NA_real_, RMSEA=NA_real_, SRMR=NA_real_,
  parameterization = parameterization_used,
  note=NA_character_,
  stringsAsFactors=FALSE
)

cfa_fit <- tryCatch({
  lavaan::cfa(
    model = model_syn,
    data  = dat_ord,
    ordered = items_use,
    estimator = "WLSMV",
    parameterization = parameterization_used,
    std.lv = TRUE
  )
}, error=function(e) e)

if(inherits(cfa_fit, "error")){
  fit_row$note <- paste0("CFA_ERROR: ", cfa_fit$message,
                         " | dropped_const=", paste(dropped_const, collapse=","))
} else {
  fm <- lavaan::fitMeasures(cfa_fit, c("chisq","df","pvalue","cfi","tli","rmsea","srmr"))
  fit_row$chisq  <- as.numeric(fm["chisq"])
  fit_row$df     <- as.numeric(fm["df"])
  fit_row$pvalue <- as.numeric(fm["pvalue"])
  fit_row$CFI    <- as.numeric(fm["cfi"])
  fit_row$TLI    <- as.numeric(fm["tli"])
  fit_row$RMSEA  <- as.numeric(fm["rmsea"])
  fit_row$SRMR   <- as.numeric(fm["srmr"])
  fit_row$note <- paste0("OK(WLSMV/", parameterization_used, ")",
                         " | dropped_const=", paste(dropped_const, collapse=","))
}

# -----------------------------
# Save plots (PDF: 2 pages)
# -----------------------------
pdf_file <- file.path(out_ega, "EGA_full_plots.pdf")
pdf(pdf_file, width=10, height=10)

qgraph::qgraph(
  ega_res$network,
  layout="spring",
  labels=items_use,
  groups=groups_list,
  color=pal[as.character(1:K)],
  legend=TRUE,
  legend.cex=0.9,
  vsize=5,
  label.cex=0.8,
  edge.width=1
)
mtext(paste0("EGA(glasso)  N=", N_use,
             "  n_items=", length(items_use),
             "  n_clusters=", K),
      side=3, line=1)

plot_rep_comm(rep_comm, pal, main=paste0("Community Replication (bootEGA)  B=", boot_iter))
dev.off()

# Save CFA fit
write.csv(fit_row, file.path(out_ega, "CFA_fit_full.csv"), row.names=FALSE)

cat("03_ega_bootega_cfa: saved PNG/PDF/CSV.\n")
