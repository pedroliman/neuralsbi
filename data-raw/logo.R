# Regenerate man/figures/logo.png from data-raw/logo-base.png.
#
# data-raw/logo-base.png is the hand-authored hex sticker artwork (hexagon
# border/fill, the density-over-network illustration) with no wordmark. It
# has no other source (it was not produced by a script in this repo), so it
# is checked in as-is and only the wordmark is regenerated here.
#
# Why this script exists: the wordmark was previously baked into the PNG by
# hand and ended up positioned low enough to clip against the hexagon's
# pointed bottom edge. The hexagon narrows to a point at the bottom, so the
# text must both clear the border vertically *and* fit the available width
# at whatever height it sits, which shrinks as the point is approached.
# safe_half_width() below encodes that constraint from measurements of the
# sticker geometry, and the script aborts rather than writing a logo that
# would clip.
#
# Run with: Rscript data-raw/logo.R (from the package root).

library(magick)
library(showtext)

base_path <- "data-raw/logo-base.png"
font_path <- "data-raw/fonts/Roboto-Bold.ttf"
out_path <- "man/figures/logo.png"

# Hexagon geometry, measured from data-raw/logo-base.png (600x692, pointy
# top/bottom, flat left/right sides). The hex is full width between the
# shoulders and narrows linearly to a point below shoulder_y.
hex_center_x <- 299.5
hex_full_width <- 517
shoulder_y <- 497
vertex_y <- 644

# Half the hexagon's interior width at row y (NA above the shoulder, where
# the sticker is already full width and this constraint doesn't bind).
safe_half_width <- function(y) {
  if (y <= shoulder_y) {
    return(hex_full_width / 2)
  }
  frac <- (vertex_y - y) / (vertex_y - shoulder_y)
  (hex_full_width / 2) * frac
}

font_add(family = "logo-bold", regular = font_path)
showtext_auto()

img <- image_read(base_path)

label <- "neuralsbi"
font_size <- 46
margin <- 30 # minimum clearance, in px, from the border on each side

draw <- image_draw(img)
## image_draw()'s y-axis increases downward (row order), the opposite of a
## normal R plot, so strheight() reports a *negative* height here; abs() it.
text_height <- abs(strheight(label, cex = font_size / 12, family = "logo-bold", font = 2))
text_width <- abs(strwidth(label, cex = font_size / 12, family = "logo-bold", font = 2))

text_bottom <- 560
text_top <- text_bottom - text_height
available_half_width <- safe_half_width(text_bottom)

if (text_width / 2 + margin > available_half_width) {
  dev.off()
  stop(
    "logo wordmark does not fit: needs ", round(text_width / 2 + margin),
    "px half-width, hexagon only offers ", round(available_half_width),
    "px at y = ", text_bottom, ". Lower font_size or raise text_bottom."
  )
}

## adj = c(0.5, 1) anchors the *top* of the string to y, so pass text_top
## rather than text_bottom.
text(
  x = hex_center_x, y = text_top,
  label, adj = c(0.5, 1),
  cex = font_size / 12, family = "logo-bold", font = 2,
  col = "white"
)
dev.off()

image_write(draw, out_path)
message("wrote ", out_path, " (wordmark half-width ", round(text_width / 2),
        "px, margin ", round(available_half_width - text_width / 2), "px)")
