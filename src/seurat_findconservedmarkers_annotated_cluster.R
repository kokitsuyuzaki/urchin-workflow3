source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Loading
load(infile)

# Check available metadata columns
cat("=== Available metadata columns ===\n")
cat("Columns:", paste(colnames(seurat.integrated@meta.data), collapse=", "), "\n\n")

# Find cluster column
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
  cat("No standard cluster column found. Using current Idents().\n")
  cat("Current Idents levels:", paste(levels(Idents(seurat.integrated)), collapse=", "), "\n\n")
} else {
  cat("Using cluster column:", cluster_column, "\n")
  Idents(seurat.integrated) <- cluster_column
}

# Diagnostic function to understand data structure
diagnose_data <- function(seurat_obj) {
  cat("=== Data Diagnosis ===\n")
  cat("Total cells:", ncol(seurat_obj), "\n")
  cat("Current identity:", length(unique(Idents(seurat_obj))), "clusters\n")
  cat("Sample groups:", length(unique(seurat_obj$sample)), "groups\n")
  
  # Show sample groups
  cat("Sample groups:", paste(unique(seurat_obj$sample), collapse=", "), "\n")
  
  # Check cluster vs sample cross-tabulation
  cat("\n=== Cluster x Sample Cross-tabulation ===\n")
  cross_tab <- table(Idents(seurat_obj), seurat_obj$sample)
  print(cross_tab)
  
  # Check which combinations have sufficient cells
  cat("\n=== Combinations with >=3 cells ===\n")
  sufficient_combinations <- which(cross_tab >= 3, arr.ind = TRUE)
  if(nrow(sufficient_combinations) > 0) {
    for(i in 1:nrow(sufficient_combinations)) {
      cluster_name <- rownames(cross_tab)[sufficient_combinations[i,1]]
      sample_name <- colnames(cross_tab)[sufficient_combinations[i,2]]
      cell_count <- cross_tab[sufficient_combinations[i,1], sufficient_combinations[i,2]]
      cat("Cluster", cluster_name, "in", sample_name, ":", cell_count, "cells\n")
    }
  }
  
  return(cross_tab)
}

# Run diagnosis
cross_tab <- diagnose_data(seurat.integrated)

# Function to check if cluster can be analyzed
check_cluster_validity <- function(cross_tab, cluster_name, min_groups = 2, min_cells = 3) {
  if(!cluster_name %in% rownames(cross_tab)) {
    return(list(valid = FALSE, reason = "Cluster not found"))
  }
  
  # Get row for this cluster
  cluster_row <- cross_tab[cluster_name, ]
  
  # Count groups with sufficient cells
  valid_groups <- names(cluster_row)[cluster_row >= min_cells]
  
  if(length(valid_groups) >= min_groups) {
    return(list(valid = TRUE, groups = valid_groups, cell_counts = cluster_row[valid_groups]))
  } else {
    return(list(valid = FALSE, 
                reason = paste("Only", length(valid_groups), "groups have >=", min_cells, "cells"),
                groups = valid_groups))
  }
}

# Get cluster information
clusters <- Idents(seurat.integrated)
unique_clusters <- levels(clusters)

# Check which clusters can be analyzed
cat("\n=== Checking cluster validity ===\n")
valid_clusters <- c()
for(cluster_name in unique_clusters) {
  check_result <- check_cluster_validity(cross_tab, cluster_name)
  
  if(check_result$valid) {
    valid_clusters <- c(valid_clusters, cluster_name)
    cat("✓ Cluster", cluster_name, " - found in", length(check_result$groups), "groups:", 
        paste(names(check_result$cell_counts), "(", check_result$cell_counts, ")", collapse=", "), "\n")
  } else {
    cat("✗ Cluster", cluster_name, " - skipping (", check_result$reason, ")\n")
    if(!is.null(check_result$groups) && length(check_result$groups) > 0) {
      cat("   Available in:", paste(check_result$groups, collapse=", "), "\n")
    }
  }
}

# Alternative approach: use FindConservedMarkers for valid clusters, FindMarkers for others
outList <- list()

if(length(valid_clusters) > 0) {
  cat("\n=== Processing valid clusters with FindConservedMarkers ===\n")
  for(cluster_name in valid_clusters) {
    cat("Processing Cluster:", cluster_name, "\n")
    
    tryCatch({
      # Try with more permissive parameters
      tmp <- FindConservedMarkers(seurat.integrated, 
                                 ident.1 = cluster_name,
                                 grouping.var = "sample",
                                 verbose = FALSE,
                                 min.pct = 0.1,          # Lower threshold
                                 logfc.threshold = 0.25,  # Lower threshold
                                 min.cells.feature = 1,   # Lower threshold
                                 min.cells.group = 1)     # Lower threshold
      
      if(nrow(tmp) > 0) {
        genenames <- rownames(tmp)
        outList[[paste0("Cluster_", cluster_name)]] <- cbind(genenames, tmp)
        cat("  - Success:", nrow(tmp), "conserved markers found\n")
      } else {
        cat("  - No conserved markers found, trying FindMarkers\n")
        
        # Fallback to FindMarkers for this cluster
        tmp2 <- FindMarkers(seurat.integrated, 
                           ident.1 = cluster_name,
                           verbose = FALSE,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
        
        if(nrow(tmp2) > 0) {
          genenames <- rownames(tmp2)
          outList[[paste0("Cluster_", cluster_name)]] <- cbind(genenames, tmp2)
          cat("  - Fallback success:", nrow(tmp2), "markers found with FindMarkers\n")
        } else {
          outList[[paste0("Cluster_", cluster_name)]] <- data.frame(genenames = character(0), stringsAsFactors = FALSE)
          cat("  - No markers found with either method\n")
        }
      }
      
    }, error = function(e) {
      cat("  - Error with FindConservedMarkers:", e$message, "\n")
      cat("  - Trying FindMarkers as fallback\n")
      
      tryCatch({
        tmp2 <- FindMarkers(seurat.integrated, 
                           ident.1 = cluster_name,
                           verbose = FALSE,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
        
        if(nrow(tmp2) > 0) {
          genenames <- rownames(tmp2)
          outList[[paste0("Cluster_", cluster_name)]] <<- cbind(genenames, tmp2)
          cat("  - Fallback success:", nrow(tmp2), "markers found\n")
        } else {
          outList[[paste0("Cluster_", cluster_name)]] <<- data.frame(genenames = character(0), stringsAsFactors = FALSE)
          cat("  - No markers found with fallback method\n")
        }
      }, error = function(e2) {
        cat("  - Fallback also failed:", e2$message, "\n")
        outList[[paste0("Cluster_", cluster_name)]] <<- data.frame(genenames = character(0), stringsAsFactors = FALSE)
      })
    })
  }
}

# Process remaining clusters with FindMarkers only
remaining_clusters <- setdiff(unique_clusters, valid_clusters)
if(length(remaining_clusters) > 0) {
  cat("\n=== Processing remaining clusters with FindMarkers ===\n")
  for(cluster_name in remaining_clusters) {
    cat("Processing Cluster:", cluster_name, "(FindMarkers only)\n")
    
    tryCatch({
      tmp <- FindMarkers(seurat.integrated, 
                        ident.1 = cluster_name,
                        verbose = FALSE,
                        min.pct = 0.1,
                        logfc.threshold = 0.25)
      
      if(nrow(tmp) > 0) {
        genenames <- rownames(tmp)
        outList[[paste0("Cluster_", cluster_name)]] <- cbind(genenames, tmp)
        cat("  - Success:", nrow(tmp), "markers found\n")
      } else {
        outList[[paste0("Cluster_", cluster_name)]] <- data.frame(genenames = character(0), stringsAsFactors = FALSE)
        cat("  - No markers found\n")
      }
      
    }, error = function(e) {
      cat("  - Error:", e$message, "\n")
      outList[[paste0("Cluster_", cluster_name)]] <- data.frame(genenames = character(0), stringsAsFactors = FALSE)
    })
  }
}

# Ensure all clusters are represented
for(cluster_name in unique_clusters) {
  sheet_name <- paste0("Cluster_", cluster_name)
  if(is.null(outList[[sheet_name]])) {
    outList[[sheet_name]] <- data.frame(genenames = character(0), stringsAsFactors = FALSE)
  }
}

# Reorder according to cluster order
expected_names <- paste0("Cluster_", unique_clusters)
outList <- outList[expected_names]

# Save
if(length(outList) > 0) {
  # Check if any results were found
  non_empty_results <- sum(sapply(outList, nrow) > 0)
  
  write_xlsx(outList, outfile)
  cat("\nResults saved to:", outfile, "\n")
  cat("Successfully processed", length(valid_clusters), "clusters with FindConservedMarkers\n")
  cat("Processed", length(remaining_clusters), "clusters with FindMarkers\n")
  cat("Total results with markers:", non_empty_results, "out of", length(unique_clusters), "clusters\n")
  
  # Show summary of results
  cat("\n=== Summary of results ===\n")
  for(sheet_name in names(outList)) {
    if(nrow(outList[[sheet_name]]) > 0) {
      if("max_pval" %in% colnames(outList[[sheet_name]])) {
        cat(sheet_name, ": ", nrow(outList[[sheet_name]]), " conserved markers\n")
      } else {
        pos_markers <- sum(outList[[sheet_name]]$avg_log2FC > 0)
        neg_markers <- sum(outList[[sheet_name]]$avg_log2FC < 0)
        cat(sheet_name, ": ", pos_markers, " up-regulated, ", neg_markers, " down-regulated markers\n")
      }
    } else {
      cat(sheet_name, ": No markers found\n")
    }
  }
  
} else {
  cat("No valid clusters found for analysis\n")
}

cat("\n=== Result interpretation ===\n")
cat("Each sheet represents one cluster's markers\n")
cat("- Conserved markers: consistent across all sample groups\n")
cat("- Regular markers: cluster vs all other clusters\n")
cat("Columns will vary depending on the method used:\n")
cat("- FindConservedMarkers: includes group-specific statistics\n")
cat("- FindMarkers: includes p_val, avg_log2FC, pct.1, pct.2, p_val_adj\n")