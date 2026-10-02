## Figure 2A, domain architecture and resistance positions.

source("R/00_theme.R")

dom    <- read.csv("data/fig1b_domains.csv",   stringsAsFactors = FALSE)
posn   <- read.csv("data/fig1b_positions.csv", stringsAsFactors = FALSE)
tracks <- read.csv("data/fig1b_tracks.csv",    stringsAsFactors = FALSE) %>%
  arrange(desc(track_order)) %>%
  mutate(y = seq_len(n()) * 1.75,
         span = to_aa - from_aa + 1)

FILLS <- c(`NS3/4A` = "amber", NS5A = "teal", NS5B = "purple")
scale_to <- function(aa, tr) (aa - tr$from_aa) / tr$span

d <- tracks %>%
  rowwise() %>%
  do({
    tr <- .
    dd <- dom %>% filter(protein == tr$protein,
                         end_aa >= tr$from_aa, start_aa <= tr$to_aa) %>%
      mutate(s = pmax(start_aa, tr$from_aa), e = pmin(end_aa, tr$to_aa),
             x1 = scale_to(s - 1, tr), x2 = scale_to(e, tr),
             y = tr$y, track_label = tr$track_label, zoom = tr$zoom,
             key = grepl("domain I$|protease|palm", domain),
             fill = ifelse(key, PAL[[FILLS[[tr$protein]]]], PAL[["grey"]]),
             line = PAL[[paste0(FILLS[[tr$protein]], "_dark")]],
             mid = (x1 + x2) / 2, wide = (x2 - x1) > 0.055)
    dd
  }) %>% ungroup()

pt <- tracks %>%
  rowwise() %>%
  do({
    tr <- .
    keep <- posn %>% filter(protein == tr$protein,
                            position >= tr$from_aa, position <= tr$to_aa)
    if (tr$protein == "NS5A" && tr$zoom == 0) keep <- keep[0, ]
    keep %>% mutate(x = scale_to(position - 0.5, tr), y = tr$y,
                    track_label = tr$track_label)
  }) %>% ungroup() %>%
  mutate(col = recode(behaviour,
                      persists = PAL[["teal_dark"]], reverts = PAL[["coral_dark"]],
                      accessory = PAL[["amber_dark"]], no_resistance = PAL[["grey_dark"]])) %>%
  arrange(y, x)

MINGAP  <- 0.046
LEAD_UP <- 0.30
LEAD_Y  <- 0.62
LAB_Y   <- 0.80

spread <- function(x, gap, lo = 0, hi = 1) {
  o <- order(x); v <- x[o]
  for (i in seq_along(v)[-1]) v[i] <- max(v[i], v[i - 1] + gap)
  over <- v[length(v)] - hi
  if (over > 0) v <- v - over
  for (i in rev(seq_along(v))[-1]) v[i] <- min(v[i], v[i + 1] - gap)
  v[v < lo] <- lo + (which(v < lo) - 1) * gap
  out <- numeric(length(x)); out[o] <- v; out
}

pt <- pt %>% group_by(y) %>% mutate(lab_x = spread(x, MINGAP)) %>% ungroup() %>%
  mutate(id = row_number())

leader <- pt %>%
  rowwise() %>%
  do(data.frame(id = .$id, col = .$col,
                lx = c(.$x, .$x, .$lab_x, .$lab_x),
                ly = c(.$y + 0.26, .$y + LEAD_UP, .$y + LEAD_Y, .$y + LEAD_Y + 0.08))) %>%
  ungroup()

zr  <- tracks %>% filter(zoom == 1)
src <- tracks %>% filter(protein == zr$protein, zoom == 0)
zx1 <- scale_to(zr$from_aa - 1, src); zx2 <- scale_to(zr$to_aa, src)
wedge <- data.frame(x = c(zx1, zx2, 1, 0),
                    y = c(src$y - 0.23, src$y - 0.23, zr$y + 0.24, zr$y + 0.24))

p <- ggplot() +
  geom_polygon(data = wedge, aes(x = x, y = y), fill = PAL[["teal"]], alpha = 0.30,
               colour = NA) +
  annotate("rect", xmin = zx1, xmax = zx2, ymin = src$y - 0.23, ymax = src$y + 0.23,
           fill = NA, colour = PAL[["teal_dark"]], linewidth = 1.2, linetype = "dashed") +
  geom_rect(data = d, aes(xmin = x1, xmax = x2 - 0.002, ymin = y - 0.22, ymax = y + 0.22,
                          fill = I(fill), colour = I(line)), linewidth = 0.85) +
  geom_text(data = filter(d, wide), aes(x = mid, y = y - 0.42, label = domain),
            size = 4.0, fontface = "bold", colour = INK) +
  geom_segment(data = pt, aes(x = x, xend = x, y = y - 0.26, yend = y + 0.26,
                              colour = I(col)), linewidth = 1.3) +
  geom_point(data = pt, aes(x = x, y = y + 0.26, colour = I(col)), size = 2.4) +
  geom_path(data = leader, aes(x = lx, y = ly, group = id, colour = I(col)),
            linewidth = 0.6, lineend = "round", linejoin = "round") +
  geom_text(data = pt, aes(x = lab_x, y = y + LAB_Y, label = position, colour = I(col)),
            size = 3.9, fontface = "bold") +
  geom_text(data = tracks, aes(x = 1.015, y = y, label = paste0(span, " aa")),
            size = 3.7, fontface = "bold", colour = INK_SOFT, hjust = 0) +
  scale_y_continuous(breaks = tracks$y, labels = tracks$track_label,
                     limits = c(min(tracks$y) - 0.75, max(tracks$y) + 1.15),
                     expand = c(0, 0)) +
  scale_x_continuous(limits = c(-0.02, 1.13), expand = c(0, 0)) +
  labs_panel("fig1b_domains",
             title = "Domain architecture and resistance positions",
             subtitle = "NS5A is shown whole and then with domain I expanded, since all ten relevant positions lie inside it",
             caption = "Teal = persists after treatment. Pink = reverts. Amber = accessory position. Grey = confers no sofosbuvir resistance.",
             x = NULL, y = NULL) +
  theme_panel() +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_text(size = BASE + 1))

save_panel(p, "fig1b_domains", w = 12, h = 7.1)
