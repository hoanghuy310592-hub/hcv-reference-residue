## Scores the reference-residue matrix against EASL 2020 Table 7 and rebuilds
## Supplementary Table S2. Run from the repository root.

suppressMessages({library(jsonlite); library(dplyr)})

POSITIONS <- c(paste0("NS5A-", c(24,28,30,31,32,38,58,62,92,93)),
               paste0("NS5B-", c(159,282,316,321)))
GENOTYPES <- paste0("GT", 1:8)

mat  <- read.csv("data/fig2_full_residue_matrix.csv", stringsAsFactors = FALSE,
                 check.names = FALSE)
easl <- fromJSON("data/parc_ras_easl2020.json", simplifyVector = FALSE)

listed <- function(pos) {
  parts <- strsplit(pos, "-")[[1]]
  genes <- if (parts[1] == "NS5A") "NS5A" else c("NS5B_nuc", "NS5B_nonnuc")
  out <- list()
  for (g in genes) {
    entry <- easl[[g]][[parts[2]]]
    if (is.null(entry)) next
    for (k in names(entry)) out[[k]] <- union(out[[k]], unlist(entry[[k]]))
  }
  out
}

gt_key <- function(gt, per_gt) {
  if (gt == "GT1") return(union(per_gt[["1a"]], per_gt[["1b"]]))
  per_gt[[substring(gt, 3)]]
}

rows <- list()
for (pos in POSITIONS) {
  per_gt <- listed(pos)
  for (gt in GENOTYPES) {
    sub <- mat[mat$genotype == gt, , drop = FALSE]
    obs <- sort(unique(sub[[pos]]))
    obs <- obs[!obs %in% c("X", "?", "-")]
    ref <- gt_key(gt, per_gt)
    if (is.null(ref) || length(ref) == 0) {
      status <- "not scored"; hit <- character(0)
    } else {
      hit <- intersect(obs, ref)
      status <- if (length(hit) > 0) "listed residue present" else "scored, none listed"
    }
    rows[[length(rows) + 1]] <- data.frame(
      position = pos, genotype = gt,
      n_subtype_references = nrow(sub),
      reference_residues_observed = paste(obs, collapse = " "),
      EASL_2020_listed_substitutions = if (is.null(ref)) "" else paste(sort(ref), collapse = " "),
      status = status,
      listed_residues_observed = paste(sort(hit), collapse = " "),
      stringsAsFactors = FALSE)
  }
}
out <- bind_rows(rows)

dir.create("data", showWarnings = FALSE)
write.csv(out, "data/Table_S2_genotype_position_matrix.csv", row.names = FALSE)

tab <- table(out$status)
evaluable <- sum(tab[c("listed residue present", "scored, none listed")])
cat("\ncells                    ", nrow(out), "\n")
cat("evaluable                ", evaluable, "\n")
cat("listed residue present   ", tab[["listed residue present"]], "\n")
cat("scored, none listed      ", tab[["scored, none listed"]], "\n")
cat("not scored               ", tab[["not scored"]], "\n")

FIG3B <- c("NS5A-30", "NS5A-31", "NS5A-62", "NS5A-93", "NS5B-282")
p <- out[out$position %in% FIG3B, ]
t2 <- table(p$status)
cat("\nFigure 3B: cells", nrow(p),
    "| evaluable", sum(t2[c("listed residue present", "scored, none listed")]),
    "| listed residue present", t2[["listed residue present"]], "\n\n")
