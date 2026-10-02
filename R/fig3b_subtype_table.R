## Figure 4B, non-epidemic subtype table.

source("R/00_theme.R")

d <- read.csv("data/fig3b_subtypes.csv", stringsAsFactors = FALSE) %>%
  mutate(row = rev(seq_len(n())),
         accent = ifelse(phenotypic_data_published == 1,
                         PAL[["teal_dark"]], PAL[["coral_dark"]]),
         band = ifelse(row %% 2 == 0, "#F7FAFB", "white"))

cols <- data.frame(
  x     = c(0.00, 0.16, 0.42, 0.68),
  head  = c("Subtype", "Distribution", "Natural polymorphism", "Reported evidence"),
  field = c("subtype", "distribution", "polymorphism", "evidence"),
  stringsAsFactors = FALSE)

lab <- do.call(rbind, lapply(seq_len(nrow(cols)), function(i)
  data.frame(x = cols$x[i] + 0.012, row = d$row,
             txt = d[[cols$field[i]]],
             bold = cols$field[i] == "subtype",
             stringsAsFactors = FALSE)))

p <- ggplot() +
  geom_rect(data = d, aes(xmin = 0, xmax = 1, ymin = row - 0.5, ymax = row + 0.5,
                          fill = I(band)), colour = NA) +
  geom_rect(data = d, aes(xmin = 0, xmax = 0.008, ymin = row - 0.45, ymax = row + 0.45,
                          fill = I(accent)), colour = NA) +
  geom_text(data = lab, aes(x = x, y = row, label = txt,
                            fontface = ifelse(bold, "bold", "plain")),
            hjust = 0, size = 3.9, colour = INK, lineheight = 0.9) +
  geom_text(data = cols, aes(x = x + 0.012, y = max(d$row) + 1, label = head),
            hjust = 0, size = 4.2, fontface = "bold", colour = "white") +
  geom_rect(aes(xmin = 0, xmax = 1, ymin = max(d$row) + 0.55, ymax = max(d$row) + 1.45),
            fill = INK, colour = NA) +
  geom_text(data = cols, aes(x = x + 0.012, y = max(d$row) + 1, label = head),
            hjust = 0, size = 4.2, fontface = "bold", colour = "white") +
  scale_x_continuous(limits = c(0, 1.02), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.3, max(d$row) + 1.55), expand = c(0, 0)) +
  labs_panel("fig3b_subtype_table", title = "Natural polymorphisms and in vitro susceptibility",
       x = NULL, y = NULL,
       caption = "Teal bar = phenotypic resistance data published. Pink bar = none located.") +
  theme_panel() +
  theme(panel.grid = element_blank(), axis.text = element_blank())

save_panel(p, "fig3b_subtype_table", w = 13, h = 5.2)
