## Palette reference sheet.

source("R/00_theme.R")

pairs <- data.frame(
  name  = c("teal", "coral", "amber", "purple", "green", "grey"),
  role  = c("default / reassuring path; NS5A; persists after treatment",
            "constrained / high-risk; reverts; fabricated; decompensated cirrhosis",
            "decision point; accessory position; NS3/4A; structural proteins",
            "NS5B",
            "trial-evidence regions (map only)",
            "neutral, unremarkable, no resistance"),
  stringsAsFactors = FALSE) %>%
  mutate(fill = unname(PAL[name]), outline = unname(PAL[paste0(name, "_dark")]),
         y = rev(seq_len(n())))

p <- ggplot(pairs) +
  geom_rect(aes(xmin = 0.02, xmax = 0.14, ymin = y - 0.34, ymax = y + 0.34,
                fill = I(fill), colour = I(outline)), linewidth = 1.4) +
  geom_text(aes(x = 0.17, y = y + 0.13, label = name),
            hjust = 0, size = 4.6, fontface = "bold", colour = INK) +
  geom_text(aes(x = 0.17, y = y - 0.17,
                label = paste0("fill ", toupper(fill), "   outline ", toupper(outline))),
            hjust = 0, size = 3.6, fontface = "bold", colour = INK_SOFT) +
  geom_text(aes(x = 0.46, y = y, label = role),
            hjust = 0, size = 3.7, colour = INK) +
  annotate("rect", xmin = 0.02, xmax = 0.14, ymin = nrow(pairs) + 0.66,
           ymax = nrow(pairs) + 1.34, fill = INK, colour = INK) +
  annotate("text", x = 0.17, y = nrow(pairs) + 1.0, hjust = 0, size = 4.4,
           fontface = "bold", colour = INK,
           label = paste0("text   ", INK, "   secondary text   ", INK_SOFT)) +
  scale_x_continuous(limits = c(0, 1.05), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.4, nrow(pairs) + 1.7), expand = c(0, 0)) +
  labs_panel("palette_swatches", title = NULL, x = NULL, y = NULL) +
  theme_panel() +
  theme(panel.grid = element_blank(), axis.text = element_blank())

save_panel(p, "palette_swatches", w = 11, h = 4.6)
