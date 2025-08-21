source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Loading
load(infile)

# Check available metadata columns
cat("=== Available metadata columns ===\n")
cat("Columns:", paste(colnames(seurat.integrated@meta.data), collapse=", "), "\n\n")

# Check current identity
cat("Current Idents():", class(Idents(seurat.integrated)), "\n")
cat("Current identity levels:", paste(levels(Idents(seurat.integrated)), collapse=", "), "\n\n")

# Set identity to seurat clusters (assuming standard Seurat workflow)
# Try different possible cluster column names
cluster_columns <- c("seurat_clusters", "RNA_snn_res.0.5", "RNA_snn_res.0.8", 
                     "clusters", "cluster", "res.0.5", "res.0.8")

cluster_column <- NULL
for(col in cluster_columns) {
  if(col %in% colnames(seurat.integrated@meta.data)) {
    cluster_column <- col
    break
  }
}

if(is.null(cluster_column)) {
  # If no standard cluster column found, use current Idents
  cat("No standard cluster column found. Using current Idents().\n")
  cat("If this is not correct, please specify the cluster column name.\n\n")
} else {
  cat("Using cluster column:", cluster_column, "\n")
  Idents(seurat.integrated) <- cluster_column
}

# Get all clusters
clusters <- Idents(seurat.integrated)
unique_clusters <- levels(clusters)

cat("=== Finding DEGs: Each cluster vs All others ===\n")
cat("Total cells:", ncol(seurat.integrated), "\n")
cat("Number of clusters:", length(unique_clusters), "\n")
cat("Clusters:", paste(unique_clusters, collapse=", "), "\n\n")

# Check cluster sizes
cluster_sizes <- table(clusters)
cat("=== Cluster sizes ===\n")
for(i in 1:length(cluster_sizes)) {
  cat("Cluster", names(cluster_sizes)[i], ":", cluster_sizes[i], "cells\n")
}
cat("\n")

# Method 1: Use FindAllMarkers (fastest, finds markers for all clusters at once)
cat("=== Using FindAllMarkers (recommended) ===\n")
tryCatch({
  all_markers <- FindAllMarkers(seurat.integrated, 
                               only.pos = FALSE,        # Include both up and down-regulated
                               min.pct = 0.1,           # Expressed in at least 10% of cells
                               logfc.threshold = 0.25,  # At least 0.25 log2 fold change
                               verbose = TRUE)
  
  cat("Success! Found", nrow(all_markers), "total markers\n")
  
  # Split by cluster
  outList <- split(all_markers, all_markers$cluster)
  
  # Add genenames column and reorder
  outList <- lapply(outList, function(df) {
    genenames <- rownames(df)
    result <- cbind(genenames, df)
    # Sort by avg_log2FC descending (most upregulated first)
    result <- result[order(result$avg_log2FC, decreasing = TRUE), ]
    return(result)
  })
  
  # Ensure all clusters are represented in the output list
  for(cluster_name in unique_clusters) {
    if(!cluster_name %in% names(outList)) {
      outList[[cluster_name]] <- data.frame(genenames = character(0))
    }
  }
  
  # Reorder according to cluster levels
  outList <- outList[unique_clusters]
  
  # Summary
  cat("\n=== Summary ===\n")
  for(cluster_name in names(outList)) {
    if(nrow(outList[[cluster_name]]) > 0) {
      pos_markers <- sum(outList[[cluster_name]]$avg_log2FC > 0)
      neg_markers <- sum(outList[[cluster_name]]$avg_log2FC < 0)
      cat("Cluster", cluster_name, ": ", pos_markers, " up-regulated, ", neg_markers, " down-regulated\n")
    } else {
      cat("Cluster", cluster_name, ": No markers found\n")
    }
  }
  
}, error = function(e) {
  cat("FindAllMarkers failed:", e$message, "\n")
  cat("Switching to individual FindMarkers approach...\n\n")
  
  # Method 2: Individual FindMarkers for each cluster
  outList <- list()
  
  for(cluster_name in unique_clusters) {
    cat("Processing Cluster:", cluster_name, "\n")
    
    tryCatch({
      # Find markers for this cluster vs all others
      markers <- FindMarkers(seurat.integrated, 
                            ident.1 = cluster_name,
                            ident.2 = NULL,  # Compare against all other cells
                            only.pos = FALSE,
                            min.pct = 0.1,
                            logfc.threshold = 0.25,
                            verbose = FALSE)
      
      if(nrow(markers) > 0) {
        genenames <- rownames(markers)
        result <- cbind(genenames, markers)
        # Sort by avg_log2FC descending
        result <- result[order(result$avg_log2FC, decreasing = TRUE), ]
        outList[[cluster_name]] <- result
        
        pos_markers <- sum(markers$avg_log2FC > 0)
        neg_markers <- sum(markers$avg_log2FC < 0)
        cat("  - Success:", pos_markers, "up-regulated,", neg_markers, "down-regulated\n")
      } else {
        outList[[cluster_name]] <- data.frame(genenames = character(0))
        cat("  - No markers found\n")
      }
      
    }, error = function(e) {
      cat("  - Error:", e$message, "\n")
      outList[[cluster_name]] <- data.frame(genenames = character(0))
    })
  }
  
  # Ensure all clusters are represented
  for(cluster_name in unique_clusters) {
    if(is.null(outList[[cluster_name]])) {
      outList[[cluster_name]] <- data.frame(genenames = character(0))
    }
  }
  
  # Reorder according to cluster levels
  outList <- outList[unique_clusters]
})

# Name the sheets properly
names(outList) <- paste0("Cluster_", names(outList))

# Save results
write_xlsx(outList, outfile)
cat("\nResults saved to:", outfile, "\n")

# Final summary
non_empty_results <- sum(sapply(outList, nrow) > 0)
total_markers <- sum(sapply(outList, nrow))
cat("Total DEGs found:", total_markers, "\n")
cat("Clusters with DEGs:", non_empty_results, "out of", length(unique_clusters), "\n")

# Show top genes for each cluster
cat("\n=== Top upregulated genes per cluster ===\n")
for(cluster_name in names(outList)) {
  if(nrow(outList[[cluster_name]]) > 0) {
    top_genes <- outList[[cluster_name]][outList[[cluster_name]]$avg_log2FC > 0, ]
    if(nrow(top_genes) > 0) {
      top_5 <- head(top_genes$genenames, 5)
      cat(cluster_name, ":", paste(top_5, collapse=", "), "\n")
    } else {
      cat(cluster_name, ": No upregulated genes\n")
    }
  } else {
    cat(cluster_name, ": No DEGs found\n")
  }
}

cat("\n=== Result interpretation ===\n")
cat("Each sheet represents one cluster vs all other clusters\n")
cat("Columns in results:\n")
cat("- genenames: Gene names\n")
cat("- p_val: Statistical significance (closer to 0 = more significant)\n")
cat("- avg_log2FC: Average log2 fold change (positive = higher in this cluster)\n")
cat("- pct.1: Percentage of cells expressing the gene in this cluster\n")
cat("- pct.2: Percentage of cells expressing the gene in other clusters\n")
cat("- p_val_adj: Adjusted p-value (multiple testing correction)\n")
cat("- cluster: Cluster ID\n")
cat("\nRows are sorted by avg_log2FC (most upregulated first)\n")
cat("\nGood cluster markers typically have:\n")
cat("- Low p_val_adj (< 0.05)\n")
cat("- High |avg_log2FC| (> 0.5 for strong markers)\n")
cat("- High pct.1 and low pct.2 for positive markers\n")