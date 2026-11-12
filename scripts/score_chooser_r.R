#!/usr/bin/env Rscript
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Score saved subsample partitions with chooseR's own code (Patterson-Cross, Levine and
# Menon, BMC Bioinformatics 2021; github.com/rbpatt2019/chooseR, GPL-3.0, at the commit
# pinned in scripts/fetch_chooser.sh). Protocol v1, scorer S3 (docs/protocol_matched_v1.md).
#
# find_matches, percent_match, group_sil and boot_median are evaluated from chooseR's
# R/pipeline.R unchanged. The tallying, scoring and silhouette steps below are copied from
# chooseR's examples/1_seurat_pipeline.R (its "Tallying", "Scoring" and "Silhouette" blocks).
# Only the clustering is ours: the partitions come from scripts/run_chooser_r.py.
#
# Usage: Rscript scripts/score_chooser_r.R <path/to/R/pipeline.R> <manifest.csv> <out.csv>
# Manifest columns: job_id, mode (tally or matrix), m (cells), b (subsamples), setting,
# labels_file (int32 full-data labels), data_file (tally: int32 m x b table with -1 for an
# absent cell, column-major; matrix: float64 m x m co-clustering frequency), dump_file
# (optional, tally only: where to write chooseR's frequency matrix for the cross-check),
# r_seed (set before each bootstrap, since the example script sets none).
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})
args <- commandArgs(trailingOnly = TRUE)
pipeline <- args[1]
manifest <- read.csv(args[2], stringsAsFactors = FALSE)
out <- args[3]

# pipeline.R begins by sourcing helper_functions.R, which needs Seurat for the clustering we
# do not use, so evaluate only the four scoring functions, as written.
keep <- c("find_matches", "percent_match", "group_sil", "boot_median")
for (e in parse(pipeline)) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && as.character(e[[2]]) %in% keep) {
    eval(e, globalenv())
  }
}
stopifnot(all(vapply(keep, exists, logical(1))))
`%>%` <- magrittr::`%>%`

status <- "ok"
msg <- ""
rows <- vector("list", nrow(manifest))
for (i in seq_len(nrow(manifest))) {
  job <- manifest[i, ]
  m <- job$m
  labels <- readBin(job$labels_file, "integer", n = m, size = 4)
  if (job$mode == "tally") {
    stopifnot(job$b == 100)  # percent_match() divides by its default n = 100
    tab <- matrix(readBin(job$data_file, "integer", n = m * job$b, size = 4), nrow = m, ncol = job$b)
    tab[tab < 0] <- NA
    results <- as.data.frame(tab)
    colnames(results) <- paste0("round_", seq_len(job$b))
    results <- tibble::add_column(tibble::as_tibble(results), cell = seq_len(m), .before = 1)
    # examples/1_seurat_pipeline.R, "Tallying"
    columns <- colnames(dplyr::select(results, -cell))
    mtchs <- matrix(0, nrow = dim(results)[1], ncol = dim(results)[1])
    for (col in columns) {
      mtchs <- Reduce("+", list(
        mtchs,
        find_matches(col, df = results)
      ))
    }
    # examples/1_seurat_pipeline.R, "Scoring"
    mtchs <- suppressWarnings(dplyr::mutate_all(
      dplyr::as_tibble(mtchs),
      function(x) dplyr::if_else(Re(x) > 0, percent_match(x), 0)
    ))
    if (!is.na(job$dump_file) && nzchar(job$dump_file)) {
      writeBin(as.vector(as.matrix(mtchs)), job$dump_file, size = 8)
    }
  } else {
    mtchs <- matrix(readBin(job$data_file, "double", n = m * m, size = 8), nrow = m, ncol = m)
  }
  # examples/1_seurat_pipeline.R, "Silhouette" and "Grouping"
  sil <- cluster::silhouette(
    x = as.numeric(as.character(unlist(labels))),
    dmatrix = (1 - as.matrix(mtchs))
  )
  grp <- group_sil(sil, job$setting)
  # the summarise step of the example script, one resolution at a time: boot_median(avg_sil)
  set.seed(job$r_seed)
  status <- "ok"
  msg <- ""
  ci <- withCallingHandlers(
    tryCatch(boot_median(grp$avg_sil), error = function(e) {
      status <<- "error"
      msg <<- conditionMessage(e)
      NULL
    }),
    warning = function(w) {
      if (status == "ok") status <<- "warning"
      msg <<- trimws(paste(msg, conditionMessage(w)))
      invokeRestart("muffleWarning")
    }
  )
  # boot::boot.ci returns NULL, with a printed note, when every bootstrap median is equal,
  # which happens when every cluster average is equal (often all 1.0 on the smallest
  # datasets). boot_median then returns NULL fields, and in the example script max(low_med)
  # becomes NA and no resolution is chosen. Recorded here as "no_interval".
  low <- if (is.null(ci) || length(ci$low_med) == 0) NA_real_ else ci$low_med
  high <- if (is.null(ci) || length(ci$high_med) == 0) NA_real_ else ci$high_med
  med <- if (is.null(ci) || length(ci$med) == 0) NA_real_ else ci$med
  if (status != "error" && is.na(low)) status <- "no_interval"
  rows[[i]] <- data.frame(job_id = job$job_id, n_clusters = nrow(grp), med = med,
                          med_direct = median(grp$avg_sil), low_med = low, high_med = high,
                          status = status, msg = msg, r_seed = job$r_seed)
}
write.csv(do.call(rbind, rows), out, row.names = FALSE)
