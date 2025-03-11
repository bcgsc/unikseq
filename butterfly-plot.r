# Load necessary libraries
library(ggplot2)
library(scales)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript script.R <data_file> <species_name>")
}

data_file <- args[1]
species_name <- args[2]

# Read data
dfa <- read.table(data_file, sep="\t", header = TRUE)

# Customize x-axis title
my_x_title <- bquote("Position of 25-mers on " * italic(.(species_name)) * " Mt genome")

# Save first plot as a PDF
pdf(paste0(species_name, "_butterfly-plot.pdf"), width=8, height=6)
ggplot(dfa, aes(fill=condition, y=value, x=position)) + 
  geom_col() + 
  labs(x=my_x_title, y="Proportion of species with 25-mers") + 
  coord_flip()
dev.off()

