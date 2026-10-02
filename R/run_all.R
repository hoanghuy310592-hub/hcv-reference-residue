## Rebuild every panel. Run from the repository root.

scripts <- c("fig1a_genome.R", "fig1b_domains.R", "fig1c_trajectories.R",
             "fig2a_residue_distribution.R", "fig2b_genotype_matrix.R",
             "fig3a_map.R", "fig3b_subtype_table.R", "fig4c_retreatment.R",
             "palette_swatches.R")
for (s in scripts) {
  message("--- ", s)
  source(file.path("R", s), echo = FALSE)
}
write_panel_text()
message("done. PDFs and 600 dpi PNGs are in out/")
