## Shared theme, palette and save helpers. Sourced by every panel script.

suppressMessages({library(ggplot2); library(dplyr); library(tidyr)})

INK      <- "#12304D"
INK_SOFT <- "#4A6076"

PAL <- c(
  teal   = "#A8D8D4", teal_dark   = "#0F6F6B",
  coral  = "#F3BFC6", coral_dark  = "#C23B50",
  amber  = "#F7DDAE", amber_dark  = "#C8871A",
  purple = "#CFC6E8", purple_dark = "#5B4B9C",
  green  = "#BFDFC6", green_dark  = "#3E7F55",
  grey   = "#DCE4EA", grey_dark   = "#7C8B99"
)

BASE <- 15

SHOW_TITLES <- FALSE
if (!exists(".PANEL_TEXT")) .PANEL_TEXT <- new.env(parent = emptyenv())
if (is.null(.PANEL_TEXT$rows)) .PANEL_TEXT$rows <- list()

labs_panel <- function(panel, title = NULL, subtitle = NULL, caption = NULL,
                       x = NULL, y = NULL) {
  .PANEL_TEXT$rows[[panel]] <- list(title = title, subtitle = subtitle, caption = caption)
  if (SHOW_TITLES) labs(title = title, subtitle = subtitle, caption = caption, x = x, y = y)
  else             labs(title = NULL,  subtitle = NULL,     caption = NULL,     x = x, y = y)
}

write_panel_text <- function(path = "out/panel_legend_text.md") {
  dir.create("out", showWarnings = FALSE)
  con <- file(path, open = "wt")
  writeLines(c("# Panel text for the figure legend", "",
               "Titles are switched off in the plots. Use these when writing the legend.", ""), con)
  for (nm in names(.PANEL_TEXT$rows)) {
    r <- .PANEL_TEXT$rows[[nm]]
    writeLines(c(paste0("## ", nm),
                 if (!is.null(r$title))    paste0("**Title.** ", r$title),
                 if (!is.null(r$subtitle)) paste0("**Subtitle.** ", r$subtitle),
                 if (!is.null(r$caption))  paste0("**Note.** ", r$caption), ""), con)
  }
  close(con); message("wrote ", path)
}

theme_panel <- function(base = BASE) {
  theme_minimal(base_size = base) +
    theme(
      text             = element_text(colour = INK, face = "bold"),
      plot.title       = element_text(size = base + 3, face = "bold", colour = INK,
                                      margin = margin(b = 6)),
      plot.subtitle    = element_text(size = base - 2, face = "plain", colour = INK_SOFT,
                                      margin = margin(b = 10)),
      plot.caption     = element_text(size = base - 3, face = "plain", colour = INK_SOFT,
                                      hjust = 0, margin = margin(t = 10)),
      axis.text        = element_text(size = base - 1, face = "bold", colour = INK),
      axis.title       = element_text(size = base, face = "bold", colour = INK),
      legend.text      = element_text(size = base - 2, face = "bold", colour = INK),
      legend.title     = element_blank(),
      legend.position  = "top",
      legend.key.size  = unit(11, "pt"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#E9EEF2", linewidth = 0.4),
      plot.margin      = margin(14, 18, 12, 14)
    )
}

save_panel <- function(p, name, w, h) {
  dir.create("out", showWarnings = FALSE)
  ggsave(file.path("out", paste0(name, ".pdf")), p, width = w, height = h,
         units = "in", device = cairo_pdf, bg = "white")
  ggsave(file.path("out", paste0(name, ".png")), p, width = w, height = h,
         units = "in", dpi = 600, bg = "white")
  message("wrote out/", name, ".pdf and .png  (", w, " x ", h, " in)")
}
