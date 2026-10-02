## Figure 4A, HCV burden, trial sites and non-epidemic subtypes.

source("R/00_theme.R")
suppressMessages(library(maps))

BINS   <- c(0, 0.25, 0.5, 1, 2, Inf)
LABELS <- c("< 0.25", "0.25 - 0.49", "0.50 - 0.99", "1.00 - 1.99", ">= 2.00")
RAMP   <- c("#E4DDF2", "#C6B8E4", "#A288CF", "#7B63B4", "#4F3F8C")
NODATA <- "#FFFFFF"

SERO_FILL <- "#EFEAE3"
SERO_KEY  <- "anti-HCV only, not comparable"

world <- map_data("world") %>% filter(region != "Antarctica")

prev_all <- read.csv("data/fig3a_prevalence.csv", stringsAsFactors = FALSE)

sero_countries <- prev_all$country[prev_all$measure == "anti-HCV"]

prev <- prev_all %>%
  filter(measure == "viraemic") %>%
  mutate(band = cut(prevalence_pct, breaks = BINS, labels = LABELS,
                    right = FALSE, include.lowest = TRUE))

trial <- read.csv("data/fig3a_trial_countries.csv", stringsAsFactors = FALSE)

PHENO   <- "phenotyped against DAAs in vitro"
NOPHENO <- "no isolate of this subtype phenotyped"

mk <- read.csv("data/fig3a_map_markers.csv", stringsAsFactors = FALSE) %>%
  mutate(grp = factor(ifelse(phenotypic_data_published == 1, PHENO, NOPHENO),
                      levels = c(PHENO, NOPHENO)))
mk$mkfill <- ifelse(mk$grp == PHENO, PAL[["teal_dark"]], PAL[["coral_dark"]])

wp <- left_join(world, prev, by = c("region" = "country"))
wp$band <- as.character(wp$band)
wp$band[is.na(wp$band)] <- "no estimate"
wp$band[wp$region %in% sero_countries] <- SERO_KEY
wp$band <- factor(wp$band, levels = c(LABELS, SERO_KEY, "no estimate"))

hatch <- NULL
if (requireNamespace("sf", quietly = TRUE) && length(sero_countries) > 0) {
  sero_poly <- world %>% filter(region %in% sero_countries)

  rings <- split(sero_poly, sero_poly$group)
  rings <- rings[vapply(rings, nrow, integer(1)) >= 4]
  geoms <- lapply(rings, function(g) {
    m <- as.matrix(g[, c("long", "lat")])
    if (!identical(m[1, ], m[nrow(m), ])) m <- rbind(m, m[1, ])
    sf::st_make_valid(sf::st_sfc(sf::st_polygon(list(m))))
  })
  shape <- sf::st_union(do.call(c, geoms))

  bb   <- as.numeric(sf::st_bbox(shape))
  STEP <- 1.5
  cs   <- seq(bb[2] - (bb[3] - bb[1]) - STEP, bb[4] + STEP, by = STEP)
  lines <- sf::st_sfc(lapply(cs, function(c0) sf::st_linestring(
    rbind(c(bb[1], bb[1] + c0 - bb[1]), c(bb[3], bb[3] + c0 - bb[1])))))

  seg <- suppressWarnings(sf::st_intersection(lines, shape))
  seg <- seg[!sf::st_is_empty(seg)]
  if (length(seg) > 0) {
    seg <- suppressWarnings(sf::st_cast(seg, "LINESTRING"))
    hatch <- do.call(rbind, lapply(seq_along(seg), function(i) {
      m <- sf::st_coordinates(seg[i])
      data.frame(x = m[-nrow(m), 1], y = m[-nrow(m), 2],
                 xend = m[-1, 1], yend = m[-1, 2])
    }))
  }
}

trial_poly <- world %>% filter(region %in% trial$country) %>%
  mutate(key = "hosted a registration trial site")

p <- ggplot() +
  geom_polygon(data = wp, aes(long, lat, group = group, fill = band),
               colour = "#C2CCD4", linewidth = 0.16) +
  {if (!is.null(hatch))
     geom_segment(data = hatch, aes(x = x, y = y, xend = xend, yend = yend),
                  colour = "#9A8F80", linewidth = 0.28, inherit.aes = FALSE)
   else NULL} +
  geom_polygon(data = trial_poly,
               aes(long, lat, group = group, colour = key),
               fill = NA, linewidth = 0.52) +
  geom_point(data = mk, aes(lon, lat, shape = grp),
             fill = mk$mkfill, size = 5.2, stroke = 1.1, colour = "white") +
  ggrepel::geom_label_repel(
    data = mk, aes(lon, lat, label = label),
    colour = mk$mkfill, size = 4.0, fontface = "bold", fill = "white",
    label.size = 0.7, box.padding = 0.55, point.padding = 0.4,
    min.segment.length = 0, segment.size = 0.5, seed = 42, max.overlaps = Inf) +

  scale_fill_manual(values = setNames(c(RAMP, SERO_FILL, NODATA),
                                      c(LABELS, SERO_KEY, "no estimate")),
                    name = "HCV viraemic prevalence, % of population",
                    drop = FALSE) +
  scale_colour_manual(values = c("hosted a registration trial site" = "#1C6B3A"),
                      name = NULL) +
  scale_shape_manual(values = setNames(c(21, 24), c(PHENO, NOPHENO)),
                     name = NULL) +
  guides(
    fill   = guide_legend(order = 1, nrow = 1, title.position = "top",
                          keywidth = unit(30, "pt"), keyheight = unit(11, "pt")),
    colour = guide_legend(order = 2, override.aes = list(fill = NA, linewidth = 1.3)),
    shape  = guide_legend(order = 3, nrow = 1, override.aes = list(
               fill = c(PAL[["teal_dark"]], PAL[["coral_dark"]]),
               colour = "white", size = 5, stroke = 1.1))) +

  coord_fixed(1.35, xlim = c(-160, 180), ylim = c(-56, 78)) +
  labs_panel("fig3a_map",
    title = "HCV burden, the trial evidence base, and the non-epidemic subtypes",
    subtitle = paste("Fill, viraemic prevalence. Hatched, seroprevalence only. Dark green outline, a country that hosted",
                     "a site in the pangenotypic registration trials. Markers, the non-epidemic subtypes of Table 3."),
    x = NULL, y = NULL,
    caption = paste(
      "Prevalence, most recent available estimate per country: 94 countries from WHO Global hepatitis report 2026",
      "(Global Health Observatory, indicator HEPATITIS_HCV_PREVALENCE_PER100, reference year 2024); 29 EU/EEA countries",
      "from Thomadakis 2023 (Lancet Reg Health Eur 36:100792, reference year 2019); USA from NHANES 2017-2020, adults.",
      "India from Polaris Observatory 2022 (viraemic, reference year 2020), because the WHO indicator carries no value for India.",
      "China, Ukraine and Iran are hatched: only anti-HCV seroprevalence was available, which counts cleared and cured infections",
      "and is not the same quantity as the viraemic bands. Countries returning <0.01% are shown as no estimate.",
      "Trial sites from the ClinicalTrials.gov records of ASTRAL-1 to -4, POLARIS-1 and -4, ENDURANCE-1 to -4,",
      "SURVEYOR-1 and -2 and EXPEDITION-1, -2 and -4; Hong Kong hosted the single Asian site of ASTRAL-1 but is too small to shade.",
      "An evidence-gap map, not a treatment-failure map: pibrentasvir retained activity against every unusual subtype tested in vitro.")) +
  theme_panel() +
  theme(panel.grid = element_blank(), axis.text = element_blank(),
        legend.position = "top", legend.box = "vertical",
        legend.box.just = "left", legend.margin = margin(b = 2),
        legend.title = element_text(size = BASE - 2, face = "bold", colour = INK))

save_panel(p, "fig3a_map", w = 13.5, h = 8.2)
