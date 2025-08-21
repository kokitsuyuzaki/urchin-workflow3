source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Loading
load(infile)

# Diagnostic function to understand data structure
diagnose_data <- function(seurat_obj) {
  cat("=== Data Diagnosis ===\n")
  cat("Total cells:", ncol(seurat_obj), "\n")
  cat("Current identity:", length(unique(Idents(seurat_obj))), "types\n")
  cat("Sample groups:", length(unique(seurat_obj$sample)), "groups\n")
  
  # Show sample groups
  cat("Sample groups:", paste(unique(seurat_obj$sample), collapse=", "), "\n")
  
  # Check celltype vs sample cross-tabulation
  cat("\n=== Celltype x Sample Cross-tabulation ===\n")
  cross_tab <- table(seurat_obj$celltype, seurat_obj$sample)
  print(cross_tab)
  
  # Check which combinations have sufficient cells
  cat("\n=== Combinations with >=3 cells ===\n")
  sufficient_combinations <- which(cross_tab >= 3, arr.ind = TRUE)
  if(nrow(sufficient_combinations) > 0) {
    for(i in 1:nrow(sufficient_combinations)) {
      celltype_name <- rownames(cross_tab)[sufficient_combinations[i,1]]
      sample_name <- colnames(cross_tab)[sufficient_combinations[i,2]]
      cell_count <- cross_tab[sufficient_combinations[i,1], sufficient_combinations[i,2]]
      cat(celltype_name, "in", sample_name, ":", cell_count, "cells\n")
    }
  }
  
  return(cross_tab)
}

# Run diagnosis
cross_tab <- diagnose_data(seurat.integrated)

# Function to check if celltype can be analyzed
check_celltype_validity <- function(cross_tab, celltype_name, min_groups = 2, min_cells = 3) {
  if(!celltype_name %in% rownames(cross_tab)) {
    return(list(valid = FALSE, reason = "Celltype not found"))
  }
  
  # Get row for this celltype
  celltype_row <- cross_tab[celltype_name, ]
  
  # Count groups with sufficient cells
  valid_groups <- names(celltype_row)[celltype_row >= min_cells]
  
  if(length(valid_groups) >= min_groups) {
    return(list(valid = TRUE, groups = valid_groups, cell_counts = celltype_row[valid_groups]))
  } else {
    return(list(valid = FALSE, 
                reason = paste("Only", length(valid_groups), "groups have >=", min_cells, "cells"),
                groups = valid_groups))
  }
}

# Set identity to celltype
Idents(seurat.integrated) <- "celltype"
celltype <- seurat.integrated@meta.data$celltype
celltype <- factor(celltype, levels=sort(unique(celltype)))

# Check which celltypes can be analyzed
cat("\n=== Checking celltype validity ===\n")
valid_celltypes <- c()
for(i in seq(nlevels(celltype))) {
  celltype_name <- levels(celltype)[i]
  check_result <- check_celltype_validity(cross_tab, celltype_name)
  
  if(check_result$valid) {
    valid_celltypes <- c(valid_celltypes, celltype_name)
    cat("✓ ", celltype_name, " - found in", length(check_result$groups), "groups:", 
        paste(names(check_result$cell_counts), "(", check_result$cell_counts, ")", collapse=", "), "\n")
  } else {
    cat("✗ ", celltype_name, " - skipping (", check_result$reason, ")\n")
    if(!is.null(check_result$groups) && length(check_result$groups) > 0) {
      cat("   Available in:", paste(check_result$groups, collapse=", "), "\n")
    }
  }
}

# Alternative approach: use FindAllMarkers for celltypes that can't use FindConservedMarkers
outList <- list()

if(length(valid_celltypes) > 0) {
  cat("\n=== Processing valid celltypes with FindConservedMarkers ===\n")
  for(celltype_name in valid_celltypes) {
    cat("Processing:", celltype_name, "\n")
    
    tryCatch({
      # Try with more permissive parameters
      tmp <- FindConservedMarkers(seurat.integrated, 
                                 ident.1 = celltype_name,
                                 grouping.var = "sample",
                                 verbose = FALSE,
                                 min.pct = 0.1,          # Lower threshold
                                 logfc.threshold = 0.25,  # Lower threshold
                                 min.cells.feature = 1,   # Lower threshold
                                 min.cells.group = 1)     # Lower threshold
      
      if(nrow(tmp) > 0) {
        genenames <- rownames(tmp)
        outList[[celltype_name]] <- cbind(genenames, tmp)
        cat("  - Success:", nrow(tmp), "conserved markers found\n")
      } else {
        cat("  - No conserved markers found, trying FindAllMarkers\n")
        
        # Fallback to FindAllMarkers for this celltype
        tmp2 <- FindMarkers(seurat.integrated, 
                           ident.1 = celltype_name,
                           verbose = FALSE,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
        
        if(nrow(tmp2) > 0) {
          genenames <- rownames(tmp2)
          outList[[celltype_name]] <- cbind(genenames, tmp2)
          cat("  - Fallback success:", nrow(tmp2), "markers found with FindMarkers\n")
        } else {
          outList[[celltype_name]] <- data.frame(genenames = character(0), stringsAsFactors = FALSE)
          cat("  - No markers found with either method\n")
        }
      }
      
    }, error = function(e) {
      cat("  - Error with FindConservedMarkers:", e$message, "\n")
      cat("  - Trying FindMarkers as fallback\n")
      
      tryCatch({
        tmp2 <- FindMarkers(seurat.integrated, 
                           ident.1 = celltype_name,
                           verbose = FALSE,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
        
        if(nrow(tmp2) > 0) {
          genenames <- rownames(tmp2)
          outList[[celltype_name]] <<- cbind(genenames, tmp2)
          cat("  - Fallback success:", nrow(tmp2), "markers found\n")
        } else {
          outList[[celltype_name]] <<- data.frame(genenames = character(0), stringsAsFactors = FALSE)
          cat("  - No markers found with fallback method\n")
        }
      }, error = function(e2) {
        cat("  - Fallback also failed:", e2$message, "\n")
        outList[[celltype_name]] <<- data.frame(genenames = character(0), stringsAsFactors = FALSE)
      })
    })
  }
}

# Process remaining celltypes with FindMarkers only
remaining_celltypes <- setdiff(levels(celltype), valid_celltypes)
if(length(remaining_celltypes) > 0) {
  cat("\n=== Processing remaining celltypes with FindMarkers ===\n")
  for(celltype_name in remaining_celltypes) {
    cat("Processing:", celltype_name, "(FindMarkers only)\n")
    
    tryCatch({
      tmp <- FindMarkers(seurat.integrated, 
                        ident.1 = celltype_name,
                        verbose = FALSE,
                        min.pct = 0.1,
                        logfc.threshold = 0.25)
      
      if(nrow(tmp) > 0) {
        genenames <- rownames(tmp)
        outList[[celltype_name]] <- cbind(genenames, tmp)
        cat("  - Success:", nrow(tmp), "markers found\n")
      } else {
        outList[[celltype_name]] <- data.frame(genenames = character(0), stringsAsFactors = FALSE)
        cat("  - No markers found\n")
      }
      
    }, error = function(e) {
      cat("  - Error:", e$message, "\n")
      outList[[celltype_name]] <- data.frame(genenames = character(0), stringsAsFactors = FALSE)
    })
  }
}

# Ensure output list has same order as original celltype levels
outList <- outList[levels(celltype)]

# Add empty dataframes for any missing celltypes
for(celltype_name in levels(celltype)) {
  if(is.null(outList[[celltype_name]])) {
    outList[[celltype_name]] <- data.frame(genenames = character(0), stringsAsFactors = FALSE)
  }
}

# Save
if(length(outList) > 0) {
  # Check if any results were found
  non_empty_results <- sum(sapply(outList, nrow) > 0)
  
  write_xlsx(outList, outfile)
  cat("\nResults saved to:", outfile, "\n")
  cat("Successfully processed", length(valid_celltypes), "celltypes with FindConservedMarkers\n")
  cat("Processed", length(remaining_celltypes), "celltypes with FindMarkers\n")
  cat("Total results with markers:", non_empty_results, "out of", nlevels(celltype), "celltypes\n")
} else {
  cat("No valid celltypes found for analysis\n")
}