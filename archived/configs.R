configs <- hash()

configs[["0.1"]] <- list(
  min.cov.bases = 10,
  qvalue.cutoff = 0.01,
  methdiff.cutoff = 10,
  log2FC.cutoff = 1,
  up.flank = 1000,
  down.flank = 1000
)

configs[["0.2"]] <- list(
  min.cov.bases = 3,
  qvalue.cutoff = 0.05,
  methdiff.cutoff = 10,
  log2FC.cutoff = 1,
  up.flank = 2000,
  down.flank = 2000
)

configs[["0.3"]] <- list(
  min.cov.bases = 10,
  qvalue.cutoff = 0.01,
  methdiff.cutoff = 10,
  log2FC.cutoff = 1,
  up.flank = 2000,
  down.flank = 1000
)

configs[["0.4"]] <- list(
  min.cov.bases = 3,
  qvalue.cutoff = 0.05,
  methdiff.cutoff = 10,
  log2FC.cutoff = 1,
  up.flank = 2000,
  down.flank = 1000
)
