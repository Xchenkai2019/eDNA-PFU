# ---- Step 1: Load required packages ----
library(corrplot)
library(ggplot2)
library(ggcorrplot)

# ---- Step 2: Read data ----
# CSV should contain correlation matrix with sample/feature names
data1 <- read.csv(
  file = "./DH-W-Genus_abound_choose_0.5_EE_矩阵新.csv",
  header = TRUE,
  row.names = "T"
)

# ---- Step 3: Convert to matrix format ----
data2 <- as.matrix(data1)

# ---- Step 4: Define custom color palette ----
mypalette <- colorRampPalette(c("#5078EF", "#FFFFFF", "#F13333"))
mycolors <- mypalette(200)  # Generate 200 gradient colors

# ---- Step 5: Create and save correlation plot ----
pdf(
  file = "DH-W-Genus_abound_choose_0.5_EE_矩阵新.pdf",
  width = 10,
  height = 10,
  family = "Times"
)

corrplot(
  data2,
  order = "original",     # Keep original variable order
  method = "color",       # Color-filled squares
  col = mycolors,         # Apply gradient colors
  type = "lower",         # Show only lower triangle
  tl.col = "black",       # Text label color
  tl.srt = 45,            # Rotate text labels
  cl.pos = "n",           # Hide color legend
  diag = TRUE             # Keep diagonal elements
)

dev.off()  # Close PDF device














