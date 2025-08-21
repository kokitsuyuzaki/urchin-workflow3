source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Loading
load(infile)

# Set identity to celltype
Idents(seurat.integrated) <- "celltype"

# Get all cell types
celltype <- seurat.integrated@meta.data$celltype
celltype <- factor(celltype, levels=sort(unique(celltype)))

cat("=== Finding markers: Each celltype vs All others ===\n")
cat("Total cells:", ncol(seurat.integrated), "\n")
cat("Cell types:", nlevels(celltype), "\n")
cat("Cell types:", paste(levels(celltype), collapse=", "), "\n\n")

# Method 1: Use FindAllMarkers (fastest, finds markers for all celltypes at once)
cat("=== Using FindAllMarkers (recommended) ===\n")
tryCatch({
  all_markers <- FindAllMarkers(seurat.integrated, 
                               only.pos = FALSE,        # Include both up and down-regulated
                               min.pct = 0.1,           # Expressed in at least 10% of cells
                               logfc.threshold = 0.25,  # At least 0.25 log2 fold change
                               verbose = TRUE)
  
  cat("Success! Found", nrow(all_markers), "total markers\n")
  
  # Split by cluster/celltype
  outList <- split(all_markers, all_markers$cluster)
  
  # Add genenames column and reorder
  outList <- lapply(outList, function(df) {
    genenames <- rownames(df)
    cbind(genenames, df)
  })
  
  # Summary
  cat("\n=== Summary ===\n")
  for(celltype_name in names(outList)) {
    pos_markers <- sum(outList[[celltype_name]]$avg_log2FC > 0)
    neg_markers <- sum(outList[[celltype_name]]$avg_log2FC < 0)
    cat(celltype_name, ": ", pos_markers, " up-regulated, ", neg_markers, " down-regulated\n")
  }
  
}, error = function(e) {
  cat("FindAllMarkers failed:", e$message, "\n")
  cat("Switching to individual FindMarkers approach...\n\n")
  
  # Method 2: Individual FindMarkers for each celltype
  outList <- list()
  
  for(i in seq(nlevels(celltype))) {
    celltype_name <- levels(celltype)[i]
    cat("Processing:", celltype_name, "\n")
    
    tryCatch({
      # Find markers for this celltype vs all others
      markers <- FindMarkers(seurat.integrated, 
                            ident.1 = celltype_name,
                            ident.2 = NULL,  # Compare against all other cells
                            only.pos = FALSE,
                            min.pct = 0.1,
                            logfc.threshold = 0.25,
                            verbose = FALSE)
      
      if(nrow(markers) > 0) {
        genenames <- rownames(markers)
        outList[[celltype_name]] <- cbind(genenames, markers)
        
        pos_markers <- sum(markers$avg_log2FC > 0)
        neg_markers <- sum(markers$avg_log2FC < 0)
        cat("  - Success:", pos_markers, "up-regulated,", neg_markers, "down-regulated\n")
      } else {
        outList[[celltype_name]] <- data.frame(genenames = character(0))
        cat("  - No markers found\n")
      }
      
    }, error = function(e) {
      cat("  - Error:", e$message, "\n")
      outList[[celltype_name]] <<- data.frame(genenames = character(0))
    })
  }
})

# Ensure all celltypes are represented
for(celltype_name in levels(celltype)) {
  if(is.null(outList[[celltype_name]])) {
    outList[[celltype_name]] <- data.frame(genenames = character(0))
  }
}

# Reorder according to original celltype levels
outList <- outList[levels(celltype)]

# Save results
write_xlsx(outList, outfile)
cat("\nResults saved to:", outfile, "\n")

# Final summary
non_empty_results <- sum(sapply(outList, nrow) > 0)
total_markers <- sum(sapply(outList, nrow))
cat("Total markers found:", total_markers, "\n")
cat("Cell types with markers:", non_empty_results, "out of", nlevels(celltype), "\n")

cat("\n=== Result interpretation ===\n")
cat("Columns in results:\n")
cat("- genenames: Gene names\n")
cat("- p_val: Statistical significance (closer to 0 = more significant)\n")
cat("- avg_log2FC: Average log2 fold change (positive = higher in this celltype)\n")
cat("- pct.1: Percentage of cells expressing the gene in this celltype\n")
cat("- pct.2: Percentage of cells expressing the gene in other celltypes\n")
cat("- p_val_adj: Adjusted p-value (multiple testing correction)\n")
cat("\nGood markers typically have:\n")
cat("- Low p_val_adj (< 0.05)\n")
cat("- High |avg_log2FC| (> 0.5 for strong markers)\n")
cat("- High pct.1 and low pct.2 for positive markers\n")
cat("- Low pct.1 and high pct.2 for negative markers\n")