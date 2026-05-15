# Script to generate UKBAnalytica hex sticker logo
# Run this script from the package root directory

library(hexSticker)
library(ggplot2)

# Create survival curve style plot for the logo
set.seed(42)
time <- seq(0, 10, by = 0.5)
surv <- exp(-0.15 * time)
df <- data.frame(time = time, surv = surv)

# Create minimal survival plot
p <- ggplot(df, aes(x = time, y = surv)) +
  geom_line(color = "#3182BD", linewidth = 1.5) +
  geom_point(data = df[c(1, 5, 10, 15, 21), ], 
             aes(x = time, y = surv), 
             color = "#E34A33", size = 2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", 
             color = "#FBB4AE", linewidth = 0.8) +
  theme_void() +
  theme_transparent()

logo_args <- list(
  subplot = p,
  package = "UKBAnalytica",
  p_size = 18,
  p_color = "#2C3E50",
  p_fontface = "bold",
  p_y = 1.4,
  s_x = 1,
  s_y = 0.75,
  s_width = 1.3,
  s_height = 0.9,
  h_fill = "#F7F7F7",
  h_color = "#3182BD",
  h_size = 1.5,
  white_around_sticker = FALSE,
  dpi = 600
)

# Generate transparent-background raster and vector versions.
do.call(sticker, c(logo_args, filename = "man/figures/logo.png"))
#do.call(sticker, c(logo_args, filename = "man/figures/logo.svg"))
#do.call(sticker, c(logo_args, filename = "man/figures/ukbanalytica-logo.svg"))
#do.call(sticker, c(logo_args, filename = "man/figures/ukbanalytica-logo.pdf"))

message("Logo saved to man/figures/logo.png, man/figures/logo.svg, and man/figures/ukbanalytica-logo.svg")
