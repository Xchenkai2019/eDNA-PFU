# ---- Load Required Packages ----
library(vegan)
library(permute)
library(lattice)
library(ape)

library(ggplot2)       # Data visualization
library(factoextra)    # PCA visualization
library(export)        # Export plots to PowerPoint
library(ggbiplot)      # Alternative PCA biplot
library(plyr)          # Data manipulation
library(scales)        # Scaling tools for ggplot2
library(grid)          # Low-level grid functions for ggplot2


otu <- read.csv("Genus_abund-DHCJ.csv",row.names = 1)
otu <- data.frame(t(otu))#读取otu数据文件
bray_dis <- vegdist(otu, method = 'bray')#根据物种组成计算样方距离，结果以 dist 数据类型存储,需要用到的函数vegdist
write.csv(as.matrix(bray_dis),"Genus_abund-DHCJ_矩阵.csv")

# ---- Step 1: If using specific variables from another file ----
chimie.new <- bray_dis

# Select and log-transform selected variables
env_log <- data.frame(
  NH3_H = log(chimie.new$NH3.H),
  TN    = log(chimie.new$TN),
  SPC   = log(chimie.new$SPC),
  CODMn = log(chimie.new$CODMn),
  Chla  = log(chimie.new$chla),
  TP    = log(chimie.new$TP),
  DO    = chimie.new$DO,
  pH    = chimie.new$pH,
  SD    = chimie.new$SD
)

# ---- Step 2: PCA ----
env_pr <- prcomp(env_log, center = TRUE, scale. = TRUE)  # PCA with centering and scaling

# Summary of PCA results
summary(env_pr)

# ---- Step 3: Scree Plot (Eigenvalues) ----
fviz_eig(env_pr, addlabels = TRUE)  # Show % variance explained
graph2ppt(file = "EW-ScreePlot.pptx", aspectr = 1)

# ---- Step 4: PCA Biplot with Grouping ----
group <- read.csv("group.csv", header = TRUE, row.names = 1)

fviz_pca_biplot(
  env_pr,
  habillage = group$group,  # Grouping variable
  palette = c("#00AFBB", "#E7B800", "#FC4E07"),
  addEllipses = TRUE,       # Add confidence ellipses
  geom = "point",
  repel = TRUE,             # Avoid overlapping labels
  col.var = "#0D0D0D"       # Arrow color
)
graph2ppt(file = "EW-PCA-Biplot.pptx", aspectr = 1)

# ---- Step 5: Custom Scree Plot with Base R ----
sdev <- env_pr$sdev
plot(
  sdev, type = "b", lty = 1, pch = 4, lwd = 1.5,
  xlab = "Principal Component Number",
  ylab = "Eigenvalue"
)
abline(h = 1, lty = 2, lwd = 1.5)  # Kaiser criterion line
graph2ppt(file = "EW-ScreePlot2.pptx", aspectr = 1)

# ---- Step 6: Calculate r-values (Loadings) ----
r <- env_pr$rotation[, 1:2]  # PC1 and PC2 loadings
PCxis <- r[, 1] + r[, 2]     # Sum of loadings for PC1 and PC2

# ---- Step 7: Plot r-values ----
data <- data.frame(PCxis = PCxis, Variable = row.names(r))
ggplot(data, aes(x = Variable, y = PCxis)) +
  geom_point(size = 5) +
  theme_bw() +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = 2, lwd = 1.5)

graph2ppt(file = "EW-r-values.pptx", aspectr = 0.5)

# ---- Step 8: Alternative PCA Plot using ggbiplot ----
# Note: Ensure 'df' contains your PCA variables and 'env$simple' contains group info
env_pr2 <- prcomp(env[, 1:11], scale = TRUE)
ggbiplot(
  env_pr2, obs.scale = 1, var.scale = 1,
  groups = env$simple,
  ellipse = TRUE
)
