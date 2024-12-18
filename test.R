# Install and load required packages
if (!requireNamespace("fgsea", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("fgsea")
}


library(fgsea)
library(ggplot2)

# Simulate ranked gene list
set.seed(42)  # For reproducibility
genes <- paste0("Gene", 1:1000)  # Simulated gene names
scores <- rnorm(1000)  # Random ranking scores (logFC, correlations, etc.)
names(scores) <- genes  # Assign gene names to scores
scores <- sort(scores, decreasing = TRUE)  # Sort by descending rank

# Simulate a gene set (e.g., KEGG B-Cell Receptor Signaling Pathway)
gene_set <- sample(genes, size = 50)  # Select 50 genes from the ranked list
numbers <- sort(as.numeric(gsub("[^0-9]", "", gene_set)))

plot(1, type = "n", xlim = c(0, max(numbers) + 2), ylim = c(-1, 1), xlab = "X", ylab = "Y")

# Add vertical lines at each number with height from 1 to -1
for (num in numbers) {
  segments(num, 1, num, -1, col = "blue", lwd = 1)  # Draw a line from (num, 1) to (num, -1)
}


# Run fgsea to calculate enrichment
fgsea_results <- fgsea(pathways = list(B_Cell_Receptor = gene_set),
                       stats = scores,
                       minSize = 15,  # Minimum size of gene set
                       maxSize = 500,  # Maximum size of gene set
                       nperm = 1000)  # Number of permutations

# View results
print(fgsea_results)

data_fsea<-plotEnrichmentData(gene_set, scores) 
# Extract the running enrichment score and plot it
plot_data <- plotEnrichment(gene_set, scores) +
  ggtitle("KEGG B-Cell Receptor Signaling Pathway") +
  theme_minimal()

print(plot_data)
