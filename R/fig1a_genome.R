## Figure 1 genome axis.

source("R/00_theme.R")

d <- read.csv("data/fig1a_genome_segments.csv", stringsAsFactors = FALSE) %>%
  arrange(order) %>%
  mutate(xend = cumsum(length_nt), xstart = xend - length_nt,
         mid = (xstart + xend) / 2,
         fill = recode(class,
                       UTR = PAL[["grey"]], structural = PAL[["amber"]],
                       other = PAL[["grey"]], target_NS3 = PAL[["amber"]],
                       target_NS5A = PAL[["teal"]], target_NS5B = PAL[["purple"]]),
         line = recode(class,
                       UTR = PAL[["grey_dark"]], structural = PAL[["amber_dark"]],
                       other = PAL[["grey_dark"]], target_NS3 = PAL[["amber_dark"]],
                       target_NS5A = PAL[["teal_dark"]], target_NS5B = PAL[["purple_dark"]]),
         wide = length_nt >= 500)

tgt <- d %>% filter(grepl("^target_", class)) %>%
  mutate(drug = recode(class,
    target_NS3  = "NS3/4A protease inhibitors\nglecaprevir, voxilaprevir",
    target_NS5A = "NS5A inhibitors\nvelpatasvir, pibrentasvir",
    target_NS5B = "NS5B nucleotide analogues\nsofosbuvir, bemnifosbuvir"),
    ylab = c(target_NS3 = -1.15, target_NS5A = -2.05, target_NS5B = -1.15)[class],
    yseg = ylab + 0.42)

p <- ggplot() +
  geom_rect(data = d, aes(xmin = xstart, xmax = xend - 12, ymin = 0, ymax = 1,
                          fill = I(fill), colour = I(line)), linewidth = 0.8) +
  geom_text(data = filter(d, wide), aes(x = mid, y = 0.5, label = segment),
            size = 4.6, fontface = "bold", colour = INK) +
  geom_text(data = filter(d, !wide), aes(x = mid, y = 1.35, label = segment),
            size = 3.5, fontface = "bold", colour = INK_SOFT, angle = 90, hjust = 0) +
  geom_segment(data = tgt, aes(x = mid, xend = mid, y = -0.12, yend = yseg),
               colour = INK_SOFT, linewidth = 0.7) +
  geom_label(data = tgt, aes(x = mid, y = ylab, label = drug, fill = I(fill)),
             size = 3.9, fontface = "bold", colour = INK, label.size = 0.7,
             label.r = unit(4, "pt"), lineheight = 0.95) +
  scale_y_continuous(limits = c(-2.7, 2.6), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0, 9600, 2000),
                     labels = scales::comma, expand = expansion(mult = c(0.02, 0.02))) +
  labs_panel("fig1a_genome", title = "Genome organisation and the three drug targets",
       subtitle = "~9.6 kb positive-sense RNA; one open reading frame encoding a ~3,000 aa polyprotein cleaved into ten proteins",
       x = "nucleotide position", y = NULL,
       caption = "Segment lengths are UNVERIFIED in this draft; check against the H77 reference (GenBank AF009606) before publication.") +
  theme_panel() +
  theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank())

save_panel(p, "fig1a_genome", w = 12, h = 5.6)
