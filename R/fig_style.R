################################################################################
## fig_style.R -- figure style shared by every figure of the paper.
##
## Sourced by R/07_figures_1_2.R and R/08_figures_3_5.R.
##   width     160 mm, the text width of an A4 page with 2.5 cm margins
##   type      Arial; every label, tick label and legend 7 pt; panel letters 8 pt bold lowercase.
##             par(mfrow) and layout() shrink text to 0.83 (2 x 2) or 0.66 (three or more rows
##             or columns), so each figure resets par(cex = 1) after setting its panels
##   lines     axes 0.45 pt, data lines about 0.8 pt; no grid, no box, L-shaped axes
##   colour    Okabe-Ito palette, distinguishable under the common colour-vision deficiencies
##   output    PNG at 600 dpi and a vector PDF with the fonts embedded, both drawn by cairo
################################################################################

FIG_W  <- 160                   # mm
FONT   <- "Arial"
PT     <- 7                     # base type size, pt
LETTER <- 8 / PT                # panel letters at 8 pt
AX_LWD <- 0.6                   # 0.6 R units = 0.45 pt
MM     <- 1 / 25.4

COL <- c(CLOSED = "#000000", ELITE2 = "#E69F00", DIRECT = "#0072B2",
         BRIDGE = "#009E73", MAS = "#D55E00")

letter <- function(l) mtext(l, side = 3, adj = 0, line = 0.35, font = 2, cex = LETTER,
                            at = grconvertX(0.01, from = "nfc", to = "user"))

## Numeric axis with a true minus sign (U+2212) on negative tick labels.
axisM <- function(side, at = axTicks(side), ...) {
  axis(side, at = at, labels = sub("^-", "−", format(at, trim = TRUE)), ...)
}

## Draw a figure to PNG (600 dpi) and to PDF with the same code.
both <- function(stem, w, h, draw) {
  png(paste0(stem, ".png"), width = w, height = h, units = "mm", res = 600,
      family = FONT, pointsize = PT, type = "cairo")
  draw(); dev.off()
  cairo_pdf(paste0(stem, ".pdf"), width = w * MM, height = h * MM,
            family = FONT, pointsize = PT)
  draw(); dev.off()
}
