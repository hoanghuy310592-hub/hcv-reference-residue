## Figure 3B, five positions resolved by genotype.

source("R/00_theme.R")

d <- read.csv("data/fig2b_genotype_matrix.csv", stringsAsFactors = FALSE) %>%
  mutate(
    state    = ifelse(scored == 0, "unscored",
                ifelse(any_ras == 1, "ras", "none")),
    position = factor(position, levels = rev(unique(position[order(pos_order)]))),
    genotype = factor(genotype, levels = paste0("GT", 1:8)),
    fill     = c(ras = PAL[["coral"]],      none = PAL[["teal"]],      unscored = "#FFFFFF")[state],
    line     = c(ras = PAL[["coral_dark"]], none = PAL[["teal_dark"]], unscored = PAL[["grey_dark"]])[state],
    ltyp     = ifelse(state == "unscored", "22", "solid"),
    txt      = ifelse(state == "unscored", PAL[["grey_dark"]], INK),
    lab      = vapply(strsplit(residues, " "), function(v) {
                 sep <- if (length(v) >= 8) " " else "  "
                 if (length(v) <= 5) paste(v, collapse = sep)
                 else paste0(paste(v[1:ceiling(length(v)/2)], collapse = sep), "\n",
                             paste(v[(ceiling(length(v)/2)+1):length(v)], collapse = sep))
               }, character(1)),
    sub      = ifelse(state == "unscored", "not scored",
                ifelse(n_residues > 1, paste0(n_residues, " residues"), "")),
    cex      = ifelse(n_residues <= 4, 5.0,
                ifelse(n_residues <= 6, 4.2,
                 ifelse(n_residues <= 8, 3.7, 3.2))),
    nudge    = ifelse(n_residues <= 5, 0.06, 0.10))

hdr <- d %>% distinct(genotype, n_subtypes)
ny  <- length(levels(d$position))

p <- ggplot(d, aes(x = genotype, y = position)) +
  geom_tile(aes(fill = I(fill), colour = I(line), linetype = I(ltyp)),
            width = 0.92, height = 0.86, linewidth = 0.9) +
  geom_text(aes(label = lab, size = I(cex), colour = I(txt),
                y = as.numeric(position) + nudge),
            fontface = "bold", lineheight = 0.82) +
  geom_text(aes(label = sub, colour = I(ifelse(state == "unscored",
                                               PAL[["grey_dark"]], INK_SOFT))),
            size = 2.9, fontface = "bold", nudge_y = -0.29) +
  geom_text(data = hdr, aes(x = genotype, y = ny + 0.75,
                            label = paste0("n=", n_subtypes)),
            size = 3.4, fontface = "bold", colour = INK_SOFT, inherit.aes = FALSE) +
  scale_y_discrete(expand = expansion(add = c(0.55, 1.15))) +
  scale_x_discrete(position = "top") +
  labs_panel("fig2b_genotype_matrix",
       title = "The same positions, resolved by genotype",
       subtitle = "Each cell lists every reference residue observed within that genotype",
       x = NULL, y = NULL,
       caption = paste("Coral, the genotype carries a residue listed at that position in EASL 2020 Table 7.",
                       "Teal, scored and no listed residue. Open dashed cell, the position is not scored for",
                       "that genotype; genotypes 7 and 8 are not scored at any position.")) +
  theme_panel() +
  theme(panel.grid = element_blank(),
        axis.text.x.top = element_text(size = BASE + 1, face = "bold"))

save_panel(p, "fig2b_genotype_matrix", w = 11.5, h = 5.8)
