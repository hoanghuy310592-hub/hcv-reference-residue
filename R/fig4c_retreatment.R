## Figure 5C, retreatment by clinical scenario.

source("R/00_theme.R")

d <- read.csv("data/fig4c_retreatment.csv", stringsAsFactors = FALSE) %>%
  mutate(row = rev(seq_len(n())),
         accent = recode(flag,
                         default = PAL[["teal_dark"]], testing = PAL[["amber_dark"]],
                         constrained = PAL[["coral_dark"]], not_failure = PAL[["grey_dark"]]),
         badge = recode(flag,
                        default = PAL[["teal"]], testing = PAL[["amber"]],
                        constrained = PAL[["coral"]], not_failure = PAL[["grey"]]),
         band = ifelse(row %% 2 == 0, "#F7FAFB", "white"),
         scen = gsub("(.{34,}?), ", "\\1\n", scenario))

cols <- data.frame(x = c(0.00, 0.31, 0.60), head = c("Scenario", "Preferred retreatment",
                                                     "Where voxilaprevir is unavailable"),
                   field = c("scen", "preferred", "fallback"), stringsAsFactors = FALSE)
lab <- do.call(rbind, lapply(seq_len(nrow(cols)), function(i)
  data.frame(x = cols$x[i] + 0.012, row = d$row, txt = d[[cols$field[i]]],
             bold = cols$field[i] == "scen", stringsAsFactors = FALSE)))

p <- ggplot() +
  geom_rect(data = d, aes(xmin = 0, xmax = 1, ymin = row - 0.5, ymax = row + 0.5,
                          fill = I(band)), colour = NA) +
  geom_rect(data = d, aes(xmin = 0, xmax = 0.008, ymin = row - 0.45, ymax = row + 0.45,
                          fill = I(accent)), colour = NA) +
  geom_text(data = lab, aes(x = x, y = row, label = txt,
                            fontface = ifelse(bold, "bold", "plain")),
            hjust = 0, size = 3.8, colour = INK, lineheight = 0.9) +
  geom_label(data = d, aes(x = 0.885, y = row, label = testing, fill = I(badge)),
             size = 3.3, fontface = "bold", colour = INK, label.size = 0.5,
             label.r = unit(7, "pt")) +
  geom_rect(aes(xmin = 0, xmax = 1, ymin = max(d$row) + 0.55, ymax = max(d$row) + 1.45),
            fill = INK, colour = NA) +
  geom_text(data = cols, aes(x = x + 0.012, y = max(d$row) + 1, label = head),
            hjust = 0, size = 4.1, fontface = "bold", colour = "white") +
  annotate("text", x = 0.885, y = max(d$row) + 1, label = "Testing",
           size = 4.1, fontface = "bold", colour = "white") +
  scale_x_continuous(limits = c(0, 1.0), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.3, max(d$row) + 1.55), expand = c(0, 0)) +
  labs_panel("fig4c_retreatment", title = "Retreatment by clinical scenario",
       x = NULL, y = NULL,
       caption = "Pink bar = decompensated cirrhosis, the only row with no protease-inhibitor option and no alternative fallback.") +
  theme_panel() +
  theme(panel.grid = element_blank(), axis.text = element_blank())

save_panel(p, "fig4c_retreatment", w = 13.5, h = 5.6)
