# ---- Step 1: Load required packages ----
library(Mfuzz)

# ---- Step 2: Load and prepare data ----
# Read CSV file with OTUs/Genus as rows and samples as columns
df <- read.csv("Genusabound_CJ5_BAC_choose.csv", row.names = 1)
df1 <- as.matrix(df)  # Convert to matrix format

# ---- Step 3: Create and standardize ExpressionSet ----
set.seed(123)  # Ensure reproducibility
mfuzz_class <- new('ExpressionSet', exprs = df1)  # Create ExpressionSet object
mfuzz_class <- standardise(mfuzz_class)           # Standardize data (mean=0, sd=1)

# ---- Step 4: Estimate fuzzification parameter m ----
m <- mestimate(mfuzz_class)
cat("Estimated m:", m, "\n")

# ---- Step 5: Evaluate cluster number using Dmin ----
# Check minimal centroid distance for different cluster numbers (2–10)
Dmin(mfuzz_class, m = m, crange = 2:10, repeats = 3, visu = TRUE)

# ---- Step 6: Perform Mfuzz clustering ----
cluster_num <- 4  # Predefined number of clusters
mfuzz_cluster <- mfuzz(mfuzz_class, c = cluster_num, m = m)

# ---- Step 7: Extract membership matrix ----
membership_matrix <- mfuzz_cluster$membership  # Membership degrees for each OTU

# ---- Step 8: Plot histogram of maximum membership values ----
max_membership <- apply(membership_matrix, 1, max)  # Max membership per OTU
hist(
  max_membership,
  breaks = 20,
  main = "Membership Strength Distribution",
  xlab = "Max Membership"
)

# ---- Step 9: Plot clusters ----
mfuzz.plot2(
  mfuzz_class,
  cl = mfuzz_cluster,
  mfrow = c(1, cluster_num),
  time.labels = colnames(df1)
)

# ---- Step 10: Build output table ----
# Group OTUs based on membership threshold (>0.5 considered 'high')
group <- ifelse(max_membership > 0.5, "h", "b")

# Extract cluster assignment
protein_cluster <- mfuzz_cluster$cluster
protein_cluster <- cbind(df[names(protein_cluster), ], cluster = protein_cluster)

# Add group classification
protein_cluster$group <- group

# ---- Step 11: Export results to CSV ----
write.csv(protein_cluster, "Genusabound_CJ5_BAC_choose_Cluster4-1.csv")