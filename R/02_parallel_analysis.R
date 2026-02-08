if(!requireNamespace("psych", quietly=TRUE)) stop("Missing package: psych")
library(psych)

dat <- read.csv(data_file)

run_pa <- function(vars, out_png){
  items <- dat[, vars, drop=FALSE]
  items <- numify_df(items)
  items <- na.omit(items)

  # Polychoric PA (uses correlation matrix)
  R <- psych::polychoric(items)$rho
  png(out_png, width=1600, height=1200, res=200)
  pa_poly <- psych::fa.parallel(R, n.obs=nrow(items), fa="fa", fm="minres", n.iter=500, plot=TRUE)
  dev.off()

  # Raw PA (no plot)
  pa_raw <- psych::fa.parallel(items, fa="fa", fm="minres", n.iter=500, plot=FALSE)

  list(
    n = nrow(items),
    poly_nfact = pa_poly$nfact,
    raw_nfact  = pa_raw$nfact
  )
}

res_full <- run_pa(
  vars_full,
  file.path(out_pa, "parallel_analysis_full.png")
)

res_pruned <- run_pa(
  vars_pruned,
  file.path(out_pa, "parallel_analysis_pruned.png")
)

out_tbl <- data.frame(
  version = c("FULL","PRUNED"),
  N = c(res_full$n, res_pruned$n),
  nfact_polychoric = c(res_full$poly_nfact, res_pruned$poly_nfact),
  nfact_raw = c(res_full$raw_nfact, res_pruned$raw_nfact)
)

write.csv(out_tbl, file.path(out_pa, "parallel_analysis_summary.csv"), row.names=FALSE)

cat("02_parallel_analysis: saved PNGs + summary CSV.\n")
