# Load required packages
library(vegan)   # For diversity index calculation
library(picante) # For phylogenetic diversity (PD) calculation

# Read OTU table with sample IDs as row names, comma-separated
# Transpose so rows are samples, columns are taxa
otu <- read.table(
  file = 'Genus_abund-DHCJ.csv',
  row.names = 1,
  sep = ',',
  stringsAsFactors = FALSE,
  check.names = FALSE,
  header = TRUE
)
otu <- t(otu)  # Transpose matrix: samples in rows, taxa in columns

# Define function to calculate alpha diversity metrics
alpha <- function(x, tree = NULL, base = exp(1)) {
  # estimateR returns richness estimates including observed richness, Chao1, ACE, etc.
  est <- estimateR(x)
  
  # Extract richness-related metrics
  Richness <- est[1, ]  # Observed species count (richness)
  Chao1 <- est[2, ]     # Chao1 richness estimator
  ACE <- est[4, ]       # ACE richness estimator
  
  # Calculate Shannon diversity index, with customizable logarithm base
  Shannon <- diversity(x, index = 'shannon', base = base)
  
  # Calculate Simpson diversity index (Gini-Simpson index = 1 - D)
  Simpson <- diversity(x, index = 'simpson')
  
  # Calculate Pielou's evenness: Shannon / log(Richness)
  # Assumes Richness > 0 to avoid log errors
  Pielou <- Shannon / log(Richness, base)
  
  # Calculate Goods coverage: 1 - (number of singletons / total counts)
  goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  # Combine all metrics into a data frame
  result <- data.frame(
    Richness = Richness,
    Shannon = Shannon,
    Simpson = Simpson,
    Pielou = Pielou,
    Chao1 = Chao1,
    ACE = ACE,
    Goods_coverage = goods_coverage
  )
  
  # If a phylogenetic tree is provided, calculate Faith's PD
  if (!is.null(tree)) {
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[, "PD"]
    names(PD_whole_tree) <- 'PD_whole_tree'
    result <- cbind(result, PD_whole_tree)
  }
  
  return(result)
}

# Calculate alpha diversity metrics on OTU table with log base 2
alpha_all <- alpha(otu, base = 2)

# Save the results to CSV without quotes
write.csv(alpha_all, 'Genus_abund-DHCJ_alpha.csv', quote = FALSE)
