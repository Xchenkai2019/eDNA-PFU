# Load required packages
library(microeco)   # Microbial ecology analysis package
library(XML)        # For handling XML data (optional, depending on needs)
library(rgexf)      # For exporting network data (optional)
library(igraph)     # Network analysis (optional)
library(export)     # Exporting plots or tables (optional)
library(ggplot2)    # Plotting package

# Read data files
otu <- read.csv("OTU-DHCJ.csv", row.names = 1, check.names = FALSE)          # OTU table, row names are sample IDs
group <- read.csv("sample_info-DHCJ.csv", row.names = 1)                      # Sample grouping information
tax <- read.csv("taxonomy-DHCJ.csv", row.names = 1)                           # Taxonomy information
color <- read.csv("OTU_color.csv")                                           # Color and taxonomy mapping

# Create microeco dataset object
dataset <- microtable$new(
  otu_table = otu,
  sample_table = group,
  tax_table = tax
)

# Print basic dataset information
print(dataset)

# Tidy the dataset for downstream analyses
dataset$tidy_dataset()

# Calculate abundance
dataset$cal_abund()

# Save abundance data to the specified folder
dataset$save_abund(dirpath = "taxa_abund1")

# Calculate alpha diversity without phylogenetic diversity (PD)
dataset$cal_alphadiv(PD = FALSE)

# Save alpha diversity results
dataset$save_alphadiv(dirpath = "alpha_diversity")

# Calculate beta diversity without using unifrac distance
dataset$cal_betadiv(unifrac = FALSE)

# Save beta diversity results
dataset$save_betadiv(dirpath = "beta_diversity")

# Set color vector; the order corresponds to rows in OTU_color.csv
mycol <- color$color

# Create transformed abundance object; taxrank should match the taxonomy level, e.g. "Phylum" or "Phylum1"
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum1", ntaxa = 20)

# Reset taxon names to display, corresponding to the phylum column in the color file
t1[["data_taxanames"]] <- color$phylum

# Set factor levels for classification to control grouping order in the plot
t1$data_abund$Classfiy <- factor(t1$data_abund$Classfiy, levels = c("Water", "PFU"))

# Plot stacked bar chart
t1$plot_bar(
  others_color = "grey70",             # Color for 'Others' category
  facet = "Site",                     # Primary facet variable
  facet2 = "Classfiy",                # Secondary facet variable
  color_values = mycol,               # Color vector
  xtext_keep = FALSE,                 # Do not show x-axis text labels
  legend_text_italic = FALSE,         # Legend text not italicized
  barwidth = 1,                      # Width of bars
  order_facet = c("East lake", "Yangtze river")  # Order of facets
)
