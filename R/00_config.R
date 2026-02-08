set.seed(123)
options(stringsAsFactors = FALSE)

# ---- Files ----
data_file   <- file.path("data", "compost_items.csv")
syntax_xlsx <- file.path("data", "cfa_model_syntax_compost.xlsx")
sheet_name  <- "16_cfa_model_syntax"

# ---- Output dirs ----
out_poly <- file.path("outputs", "01_polychoric")
out_pa   <- file.path("outputs", "02_parallel_analysis")
out_ega  <- file.path("outputs", "03_ega")
out_cfa  <- file.path("outputs", "04_cfa_6models")
out_sum  <- file.path("outputs", "05_cfa_metrics")

dirs <- c(out_poly, out_pa, out_ega, out_cfa, out_sum)
for(d in dirs) if(!dir.exists(d)) dir.create(d, recursive = TRUE)

# ---- Items ----
vars_full <- c(
  "int1","int2","int3",
  "att1","att2","att3","att4","att5","att6",
  "sn1","sn2","sn3","sn4","sn5","sn6","sn7","sn8","sn9","sn10","sn11",
  "pbc1","pbc2","pbc3","pbc4","pbc5","pbc6","pbc7","pbc8","pbc9","pbc10","pbc11"
)

# Pruned example (remove att5, att6, pbc1)
vars_pruned <- c(
  "int1","int2","int3",
  "att1","att2","att3","att4",
  "sn1","sn2","sn3","sn4","sn5","sn6","sn7","sn8","sn9","sn10","sn11",
  "pbc2","pbc3","pbc4","pbc5","pbc6","pbc7","pbc8","pbc9","pbc10","pbc11"
)

# ---- Analysis settings ----
boot_iter <- 500
parameterization_used <- "delta"  # "delta" or "theta"

# ---- 6-model names (must exist in Excel 'model' column) ----
need_models <- c("FULL_TPB","FULL_EFA","FULL_EGA","PRUNED_TPB","PRUNED_EFA","PRUNED_EGA")
