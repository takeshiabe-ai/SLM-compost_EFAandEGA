pkgs <- c("psych","ggcorrplot","ggplot2")
invisible(lapply(pkgs, function(p){
  if(!requireNamespace(p, quietly=TRUE)) stop("Missing package: ", p)
}))
library(psych); library(ggcorrplot); library(ggplot2)

dat <- read.csv(data_file)
items <- dat[, vars_full, drop=FALSE]
items <- numify_df(items)
items <- na.omit(items)

R_poly <- psych::polychoric(items)$rho

png(file.path(out_poly, "polychoric_heatmap_full.png"), width=1200, height=900)
p <- ggcorrplot(R_poly, hc.order=FALSE, lab=FALSE, type="lower",
                outline.col="white", title="Polychoric correlations (FULL)")
print(p + ggplot2::theme_minimal(base_size=18) +
        ggplot2::theme(plot.title=element_text(face="bold", hjust=0.5),
                       axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)))
dev.off()

write.csv(R_poly, file.path(out_poly, "polychoric_full.csv"), row.names=TRUE)

cat("01_polychoric: saved heatmap + matrix.\n")
