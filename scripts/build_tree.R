#!/usr/bin/env Rscript
# =============================================================================
#  BioRemmer v1.0 — build_tree.R
#  Builds a Neighbor-Joining phylogenomic tree from a VAMPhyRE distance matrix
#  (.meg format) using ape::nj().
#  Replaces MEGA 11 (megacc) — no manual installation required.
#
#  Usage:
#    Rscript build_tree.R <input.meg> <output.nwk>
#
#  Dependencies (all in conda bioremmer environment):
#    - r-ape (installed as ggtree dependency)
# =============================================================================

suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(phangorn))  # needed for midpoint()

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  cat("Usage: Rscript build_tree.R <input.meg> <output.nwk>\n")
  quit(status = 1)
}

meg_file <- args[1]
nwk_file <- args[2]

if (!file.exists(meg_file)) {
  stop("Input file not found: ", meg_file)
}

# ── Parse MEGA distance matrix (.meg format) ─────────────────────────────────
# MEGA lower-triangular format:
#   #mega
#   !Title ...;
#   !Format DataType=Distance ...;
#   #taxon1
#   #taxon2
#   0.123
#   #taxon3
#   0.456  0.789
#   ...

parse_mega_dist <- function(filepath) {
  lines <- readLines(filepath)

  # Extract taxon names (lines starting with #, skip header #mega)
  taxon_lines <- grep("^#", lines, value = TRUE)
  # First line is "#mega" header — remove it
  taxon_names <- taxon_lines[!grepl("^#mega", taxon_lines, ignore.case = TRUE)]
  taxon_names <- sub("^#", "", taxon_names)   # strip leading #
  taxon_names <- trimws(taxon_names)
  n <- length(taxon_names)

  if (n < 3) stop("Need at least 3 taxa to build a tree. Found: ", n)

  # Extract numeric data lines (lines with numbers only)
  data_lines <- lines[grepl("^[0-9\\. \\t-]+$", trimws(lines))]
  data_lines <- data_lines[trimws(data_lines) != ""]

  # Build lower-triangular distance matrix
  dist_mat <- matrix(0, nrow = n, ncol = n,
                     dimnames = list(taxon_names, taxon_names))

  row_idx <- 2   # start filling from row 2 (row 1 has no values)
  for (i in seq_along(data_lines)) {
    vals <- as.numeric(strsplit(trimws(data_lines[[i]]), "\\s+")[[1]])
    for (j in seq_along(vals)) {
      dist_mat[row_idx, j] <- vals[j]
      dist_mat[j, row_idx] <- vals[j]   # mirror
    }
    row_idx <- row_idx + 1
    if (row_idx > n) break
  }

  as.dist(dist_mat)
}

# ── Build NJ tree ─────────────────────────────────────────────────────────────
cat(sprintf("[build_tree.R] Reading: %s\n", meg_file))
dist_matrix <- tryCatch(
  parse_mega_dist(meg_file),
  error = function(e) stop("Failed to parse .meg file: ", e$message)
)

cat(sprintf("[build_tree.R] Building NJ tree for %d taxa...\n",
            attr(dist_matrix, "Size")))

tree <- nj(dist_matrix)

# Root the tree at midpoint for cleaner visualization
tree <- midpoint(tree)

# ── Write Newick output ───────────────────────────────────────────────────────
write.tree(tree, file = nwk_file)
cat(sprintf("[build_tree.R] Tree written to: %s\n", nwk_file))
cat(sprintf("[build_tree.R] Tips: %d | Internal nodes: %d\n",
            length(tree$tip.label), tree$Nnode))
