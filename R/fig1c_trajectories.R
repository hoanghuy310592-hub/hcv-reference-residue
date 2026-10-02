## Figure 2B, three trajectories of a variant.

source("R/00_theme.R")

d <- read.csv("data/fig1c_trajectories.csv", stringsAsFactors = FALSE) %>%
  mutate(trajectory = factor(trajectory,
    levels = c("selected, persisting", "drug-relevant, costly", "lineage polymorphism")))

cols <- c("selected, persisting"  = PAL[["teal_dark"]],
          "drug-relevant, costly" = PAL[["coral_dark"]],
          "lineage polymorphism"  = PAL[["grey_dark"]])

X0     <- -0.10
X1     <-  1.00
GUTTER <-  1.62
LAB_X  <-  1.115
TIP    <-  1.045

lab <- data.frame(
  yend = c(97, 92, 3),
  ylab = c(102, 80, 6),
  text = c("lineage\npolymorphism",
           "selected,\npersisting",
           "drug-relevant,\ncostly"),
  stringsAsFactors = FALSE
)
lab$col <- c(PAL[["grey_dark"]], PAL[["teal_dark"]], PAL[["coral_dark"]])

THR_Y <- 31

p <- ggplot(d, aes(x = x, y = frequency, colour = trajectory, group = trajectory)) +

  annotate("rect", xmin = X0, xmax = 1.06, ymin = 0, ymax = 15,
           fill = PAL[["grey"]], alpha = 0.55) +
  annotate("segment", x = X0, xend = 1.06, y = 15, yend = 15,
           linetype = "dashed", colour = INK_SOFT, linewidth = 0.8) +

  annotate("segment", x = 1.06, xend = TIP, y = 15, yend = 15,
           linetype = "dashed", colour = INK_SOFT, linewidth = 0.8) +
  annotate("segment", x = TIP, xend = TIP, y = 15, yend = THR_Y,
           linetype = "dashed", colour = INK_SOFT, linewidth = 0.8) +
  annotate("text", x = LAB_X, y = THR_Y, hjust = 0, vjust = 0.5,
           size = 4.0, fontface = "bold", lineheight = 0.92, colour = INK_SOFT,
           label = "15% conventional\nreporting threshold") +

  geom_segment(data = lab, inherit.aes = FALSE,
               aes(x = X1, xend = TIP, y = yend, yend = ylab),
               colour = lab$col, linewidth = 0.9, alpha = 0.8) +

  geom_line(linewidth = 2.2, lineend = "round") +
  geom_point(size = 5.2) +

  annotate("label", x = X0 + 0.03, y = 108, hjust = 0, vjust = 1,
           label = "Schematic", size = 4.3, fontface = "bold",
           colour = INK_SOFT, fill = "white", label.size = 0,
           label.padding = unit(3, "pt")) +

  geom_text(data = lab, inherit.aes = FALSE,
            aes(x = LAB_X, y = ylab, label = text),
            colour = lab$col, hjust = 0, vjust = 0.5, size = 4.5,
            fontface = "bold", lineheight = 0.92) +

  scale_colour_manual(values = cols, guide = "none") +
  scale_x_continuous(breaks = c(0, 1),
                     labels = c("before\ntreatment", "at virological\nfailure"),
                     limits = c(X0 - 0.02, GUTTER), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 112),
                     breaks = c(4, 15, 96),
                     labels = c("low", "15%\nthreshold", "high"),
                     expand = c(0, 0)) +
  labs_panel("fig1c_trajectories",
       title = "Three trajectories of a variant",
       subtitle = "Schematic. No measured data.",
       x = NULL, y = "variant frequency\n(not to scale)",
       caption = "A single post-failure sample cannot distinguish the first trajectory from the third.") +
  theme_panel() +
  theme(legend.position  = "none",
        panel.grid.major.x = element_blank(),
        axis.text.x      = element_text(size = BASE, face = "bold", lineheight = 0.95,
                                        margin = margin(t = 6)),
        axis.title.y     = element_text(margin = margin(r = 10)),
        axis.text.y      = element_text(size = BASE - 1, face = "bold",
                                        lineheight = 0.92, colour = INK),
        panel.grid.major.y = element_blank(),
        plot.margin      = margin(16, 14, 14, 18))

save_panel(p, "fig1c_trajectories", w = 10.2, h = 5.6)
