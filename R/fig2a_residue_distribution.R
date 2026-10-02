## Figure 3A, reference residue distribution across 93 genomes.

source("R/00_theme.R")

d <- read.csv("data/fig2a_residues.csv", stringsAsFactors = FALSE)

mat <- read.csv("data/fig2_full_residue_matrix.csv", stringsAsFactors = FALSE,
                check.names = FALSE)
h77_row <- mat[mat$subtype == "1a", , drop = FALSE]
if (nrow(h77_row) != 1L) stop("expected exactly one genotype 1a row in fig2_full_residue_matrix.csv")

h77_for <- function(pos) {
  if (!pos %in% names(h77_row)) stop("position not in the residue matrix: ", pos)
  as.character(h77_row[[pos]][1])
}

d <- d %>%
  arrange(pos_order, desc(n)) %>%
  group_by(position) %>%
  mutate(xend = cumsum(n), xstart = xend - n,
         mid = (xstart + xend) / 2,
         wide = n >= 6) %>%
  ungroup() %>%
  mutate(position = factor(position, levels = rev(unique(d$position[order(d$pos_order)]))),
         fill = ifelse(is_ras == 1, PAL[["coral"]], PAL[["teal"]]),
         line = ifelse(is_ras == 1, PAL[["coral_dark"]], PAL[["teal_dark"]]))

ndist <- d %>% filter(!residue %in% c("X", "?")) %>%
  count(position, name = "n_res")

REF_X <- -6.5
href <- data.frame(position = levels(d$position), stringsAsFactors = FALSE)
href$residue <- vapply(href$position, h77_for, character(1))
href$position <- factor(href$position, levels = levels(d$position))

p <- ggplot(d) +
  geom_rect(aes(xmin = xstart, xmax = xend - 0.25,
                ymin = as.numeric(position) - 0.40,
                ymax = as.numeric(position) + 0.40,
                fill = I(fill), colour = I(line)), linewidth = 0.7) +
  geom_text(data = filter(d, wide),
            aes(x = mid, y = as.numeric(position),
                label = paste0(residue, "  ", n)),
            size = 4.3, fontface = "bold", colour = INK) +
  geom_label(data = href,
             aes(x = REF_X, y = as.numeric(position), label = residue),
             size = 4.5, fontface = "bold", colour = INK,
             fill = PAL[["grey"]], label.size = 0, label.r = unit(2.5, "pt"),
             label.padding = unit(3.2, "pt")) +
  annotate("text", x = REF_X, y = length(levels(d$position)) + 0.9,
           label = "H77", size = 3.8, fontface = "bold", colour = INK_SOFT) +

  geom_text(data = ndist, aes(x = 96, y = as.numeric(position), label = n_res),
            size = 4.4, fontface = "bold", colour = INK, hjust = 0) +
  annotate("text", x = 96, y = length(levels(d$position)) + 0.9,
           label = "distinct", size = 3.8, fontface = "bold",
           colour = INK_SOFT, hjust = 0) +
  scale_y_continuous(breaks = seq_along(levels(d$position)),
                     labels = levels(d$position),
                     expand = expansion(add = c(0.6, 1.4))) +
  scale_x_continuous(limits = c(REF_X - 4, 108), breaks = c(0, 25, 50, 75, 93),
                     expand = expansion(add = c(0.5, 0))) +
  labs_panel("fig2a_residue_distribution", title = "Reference amino acid at each scored resistance position",
       subtitle = "93 confirmed HCV subtype reference genomes; segment width is the number of subtypes carrying that residue",
       x = "number of subtype reference genomes", y = NULL,
       caption = paste("Coral, the residue is listed as a resistance substitution at that position in EASL 2020 Table 7",
                       "for at least one genotype; the panel pools all 93 genomes, so a coral flag is not genotype-specific",
                       "(see panel B). The H77 column gives the residue at that position in the genotype 1a reference,",
                       "which is the residue the standard two-letter nomenclature treats as wild type.",
                       "Resistance lists taken from the author's extraction of EASL 2020 Table 7.")) +
  theme_panel() +
  theme(panel.grid.major.y = element_blank())

save_panel(p, "fig2a_residue_distribution", w = 11.5, h = 8.2)
