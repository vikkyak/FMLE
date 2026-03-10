suppressPackageStartupMessages({
  library(optparse)
  library(here)
})

option_list <- list(
  make_option("--out_root", type="character", default=here::here("results")),
  make_option("--data_root", type="character", default=here::here("data")),
  make_option("--seed", type="integer", default=42),

  make_option("--scale_factor", type="double", default=1e4),
  make_option("--hvg_n", type="integer", default=2000),
  make_option("--pca_npcs", type="integer", default=25),
  make_option("--gate_pcs", type="integer", default=20),

  make_option("--adt_total_q", type="double", default=0.995),
  make_option("--cap_q", type="double", default=0.995)
)

opt <- parse_args(OptionParser(option_list = option_list))

cfg <- list(
  out_root      = normalizePath(opt$out_root, mustWork = FALSE),
  data_root     = normalizePath(opt$data_root, mustWork = FALSE),
  seed          = opt$seed,
  scale_factor  = opt$scale_factor,
  hvg_n         = opt$hvg_n,
  pca_npcs      = opt$pca_npcs,
  gate_pcs      = opt$gate_pcs,
  adt_total_q   = opt$adt_total_q,
  cap_q         = opt$cap_q
)

dir.create(cfg$out_root, recursive = TRUE, showWarnings = FALSE)
set.seed(cfg$seed)

