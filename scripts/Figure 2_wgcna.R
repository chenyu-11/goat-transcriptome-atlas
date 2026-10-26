#!/usr/bin/env Rscript

# ============================================================
# Complete WGCNA Analysis Script (with resume support)
# Function: Weighted Gene Co-expression Network Analysis
# Input: tpm0.1_filtered_gene_expression.csv, trait.txt
# Output: moduleTraitCor.csv, moduleTraitPvalue.csv, gene_module.csv
# Features: Automatically saves intermediate results, supports resuming from interruptions
# ============================================================

# Load required packages
suppressMessages(library(WGCNA))
suppressMessages(library(dynamicTreeCut))
suppressMessages(library(fastcluster))

# Set random seed
set.seed(12345)

# Set parameters (modify as needed)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()  # Enable multi-threading acceleration

cat("\n========================================\n")
cat("     WGCNA Analysis Started (Resume Mode)\n")
cat("========================================\n\n")

# ============================================================
# Helper functions: check whether a step is completed
# ============================================================
check_step <- function(step_name, check_file) {
  if(file.exists(check_file)) {
    cat("✅ Step '", step_name, "' already completed, skipping\n", sep = "")
    return(TRUE)
  } else {
    cat("⏳ Step '", step_name, "' not completed, starting...\n", sep = "")
    return(FALSE)
  }
}

mark_step <- function(step_name, marker_file) {
  # Create marker file
  writeLines(paste(step_name, "completed at", date()), marker_file)
  cat("✅ Step '", step_name, "' completed\n", sep = "")
}

# ============================================================
# 1. Data loading and preprocessing
# ============================================================
cat("\n=== 1. Data loading and preprocessing ===\n")

if(check_step("Data loading and preprocessing", "step1_done.txt")) {
  # If completed, load saved data
  load("WGCNA_step1.RData")
} else {
  # Read data
  gene <- read.csv("tpm0.1_filtered_gene_expression.csv", row.names = 1)
  cat("Original matrix dimensions:", dim(gene), "\n")
  
  # Transpose matrix (samples in rows, genes in columns)
  gene <- t(gene)
  cat("Transposed matrix dimensions:", dim(gene), "\n")
  
  # log2(TPM + 1) transformation
  cat("\n=== log2(TPM+1) transformation ===\n")
  gene <- log2(gene + 1)
  cat("Expression value range after transformation:", range(gene), "\n")
  cat("Matrix dimensions after transformation:", dim(gene), "\n")
  
  # Data quality control
  cat("\n=== Data quality control ===\n")
  gsg <- goodSamplesGenes(gene, verbose = 3)
  if(!gsg$allOK) {
    cat("⚠️ Low-quality genes/samples found, removing...\n")
    gene <- gene[gsg$goodSamples, gsg$goodGenes]
    cat("Matrix dimensions after filtering:", dim(gene), "\n")
  }
  
  # Save step 1 results
  save(gene, file = "WGCNA_step1.RData")
  mark_step("Data loading and preprocessing", "step1_done.txt")
}

# ============================================================
# 2. Soft threshold selection
# ============================================================
cat("\n=== 2. Soft threshold selection ===\n")

if(check_step("Soft threshold selection", "step2_done.txt")) {
  load("WGCNA_step2.RData")
} else {
  powers <- 1:20
  sft <- pickSoftThreshold(gene, powerVector = powers, verbose = 8)
  
  # ============================================================
  # Select the lowest power satisfying the scale-free topology threshold
  # ============================================================
  r2_cutoff <- 0.85
  signed_R2 <- -sign(sft$fitIndices[,3]) * sft$fitIndices[,2]
  valid_powers <- which(signed_R2 >= r2_cutoff)
  
  if(length(valid_powers) > 0) {
    best_power <- powers[min(valid_powers)]
  } else {
    best_power <- powers[which.max(signed_R2)]
    cat("⚠️ R² > 0.85 not reached, selecting power with highest R²:", best_power, "\n")
  }
  
  cat("Final selected Power:", best_power, "\n")
  
  # Save results
  save(sft, best_power, powers, signed_R2, file = "WGCNA_step2.RData")
  mark_step("Soft threshold selection", "step2_done.txt")
}

# ============================================================
# Plot soft threshold selection figure (redrawn each run)
# ============================================================
cat("\n=== Plotting soft threshold figure ===\n")

pdf("1_SoftThreshold.pdf", width = 10, height = 7)
par(mfrow = c(1, 2))

# Scale-free topology fit
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     type = "n", xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = "Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels = powers, col = "red")
abline(h = 0.85, col = "red", lwd = 2)

# Mean connectivity
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     type = "n", xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity", main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, col = "red")

dev.off()
cat("✅ Soft threshold figure saved: 1_SoftThreshold.pdf\n")

# ============================================================
# 3. Network construction and module identification
# ============================================================
cat("\n=== 3. Network construction and module identification ===\n")

if(check_step("Network construction", "step3_done.txt")) {
  load("WGCNA_step3.RData")
} else {
  # Obtain TOM matrix
  adjacency <- adjacency(gene, power = best_power)
  tom_sim <- TOMsimilarity(adjacency)
  rownames(tom_sim) <- rownames(adjacency)
  colnames(tom_sim) <- colnames(adjacency)
  cat("TOM matrix dimensions:", dim(tom_sim), "\n")
  
  # TOM dissimilarity = 1 - TOM similarity
  tom_dis <- 1 - tom_sim
  
  # Hierarchical clustering tree
  geneTree <- hclust(as.dist(tom_dis), method = 'average')
  
  # Use dynamic tree cut to mine modules
  minModuleSize <- 30
  dynamicMods <- cutreeDynamic(dendro = geneTree, distM = tom_dis,
                               deepSplit = 2, pamRespectsDendro = FALSE, 
                               minClusterSize = minModuleSize)
  cat("Number of modules:", length(unique(dynamicMods)), "\n")
  print(table(dynamicMods))
  
  # Assign module colors
  dynamicColors <- labels2colors(dynamicMods)
  
  # Plot gene clustering tree and module colors
  pdf("2_Module_Dendrogram.pdf", width = 12, height = 8)
  plotDendroAndColors(geneTree, dynamicColors, 'Dynamic Tree Cut',
                      dendroLabels = FALSE, addGuide = TRUE, hang = 0.03, 
                      guideHang = 0.05, main = 'Gene dendrogram with module colors')
  dev.off()
  cat("✅ Module clustering figure saved: 2_Module_Dendrogram.pdf\n")
  
  # Gene expression clustering tree and co-expression topology heatmap
  pdf("3_Network_Heatmap.pdf", width = 10, height = 10)
  plot_sim <- -(1 - tom_sim)
  diag(plot_sim) <- NA
  TOMplot(plot_sim, geneTree, dynamicColors,
          main = 'Network heatmap plot')
  dev.off()
  cat("✅ Network heatmap saved: 3_Network_Heatmap.pdf\n")
  
  # Save step 3 results
  save(adjacency, tom_sim, tom_dis, geneTree, dynamicMods, dynamicColors,
       file = "WGCNA_step3.RData")
  mark_step("Network construction", "step3_done.txt")
}


# ============================================================
# 4. Module eigengene calculation and merging
# ============================================================
cat("\n=== 4. Module eigengene calculation and merging ===\n")

if(check_step("Module merging", "step4_done.txt")) {
  load("WGCNA_step4.RData")
} else {
  # Calculate module eigengenes
  MEList <- moduleEigengenes(gene, colors = dynamicColors)
  MEs <- MEList$eigengenes
  
  # Module eigengene clustering
  ME_cor <- cor(MEs)
  
  # Plot clustering tree (fix y-axis display issue)
  pdf("4_Module_Eigengene_Clustering.pdf", width = 8, height = 6)
  METree <- hclust(as.dist(1 - ME_cor), method = 'average')
  
  # Adjust margins: left, bottom, top, right (increase left margin to display y-axis labels)
  par(mar = c(5, 2, 2, 2))  # increased from original c(0.4, 0.4, 0.4, 0.4)
  
  plot(METree, 
       main = 'Clustering of module eigengenes', 
       xlab = 'Modules',
       ylab = 'Height (1 - Correlation)',
       sub = '',
       cex.lab = 1.2,
       cex.main = 1.3,
       cex.axis = 0.9)

  # Red dashed line marking cut height
  abline(h = 0.25, col = 'red', lwd = 2, lty = 2)

  # Add legend
  legend("topright", 
         legend = c("cutHeight = 0.25"),
         col = "red", 
         lty = 2, 
         lwd = 2,
         bty = "n",
         cex = 0.9)
  
    dev.off()
    cat("✅ Module eigengene clustering figure saved: 4_Module_Eigengene_Clustering.pdf\n")
  
  # Merge highly similar modules
  merge_module <- mergeCloseModules(gene, dynamicColors, cutHeight = 0.25, verbose = 3)
  mergedColors <- merge_module$colors
  module <- merge_module$newMEs
  
  cat("Number of modules before merging:", length(unique(dynamicColors)), "\n")
  cat("Number of modules after merging:", length(unique(mergedColors)), "\n")
  
  pdf("5_Merged_Module_Dendrogram.pdf", width = 12, height = 8)
  plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), 
                      c('Dynamic Tree Cut', 'Merged dynamic'),
                      dendroLabels = FALSE, addGuide = TRUE, hang = 0.03, 
                      guideHang = 0.05)
  dev.off()
  cat("✅ Merged module clustering figure saved: 5_Merged_Module_Dendrogram.pdf\n")
  
  # Save step 4 results
  save(MEList, MEs, ME_cor, METree, merge_module, mergedColors, module,
       file = "WGCNA_step4.RData")
  mark_step("Module merging", "step4_done.txt")
}


# ============================================================
# 5. Export gene-module correspondence (save first for later use)
# ============================================================
cat("\n=== 5. Saving gene-module correspondence ===\n")

gene_module <- data.frame(gene_name = colnames(gene), 
                          module = mergedColors, 
                          stringsAsFactors = FALSE)
write.csv(gene_module, "gene_module.csv", row.names = FALSE)
cat("✅ Gene-module correspondence saved: gene_module.csv\n")

# Save complete intermediate results (for later recovery)
save(gene, mergedColors, module, geneTree, tom_sim, tom_dis, 
     adjacency, dynamicColors, MEList,
     file = "WGCNA_intermediate.RData")
cat("✅ Complete intermediate results saved: WGCNA_intermediate.RData\n")

# ============================================================
# 6. Module-trait association analysis (if trait.txt exists)
# ============================================================
cat("\n=== 6. Module-trait association analysis ===\n")

if(file.exists("trait.txt")) {
  cat("✅ Found trait data: trait.txt\n")
  
  # Read trait data
  trait <- read.delim('trait.txt', row.names = 1, check.names = FALSE)
  cat("Trait data dimensions:", dim(trait), "\n")
  
  # Check whether samples match
  common_samples <- intersect(rownames(module), rownames(trait))
  
  if(length(common_samples) == 0) {
    cat("⚠️ Warning: module and trait data have no common sample names\n")
    cat("  Module samples (first 5):", paste(head(rownames(module)), collapse=", "), "\n")
    cat("  Trait samples (first 5):", paste(head(rownames(trait)), collapse=", "), "\n")
    cat("  Skipping trait association analysis\n")
    
  } else {
    cat("Number of common samples:", length(common_samples), "\n")
    
    # Match samples
    module_matched <- module[common_samples, , drop = FALSE]
    trait_matched <- trait[common_samples, , drop = FALSE]
    
    # Correlation analysis between gene co-expression modules and traits
    moduleTraitCor <- cor(module_matched, trait_matched, use = 'p')
    moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(module_matched))
    
    # ============================================================
    # Plot correlation heatmap (R value on top, P value in parentheses below)
    # ============================================================
    pdf("6_Module_Trait_Heatmap.pdf", width = 12, height = 10)
    
    # Create display matrix: R value on new line with P value in parentheses
    textMatrix <- matrix("", nrow = nrow(moduleTraitCor), ncol = ncol(moduleTraitCor))
    for(i in 1:nrow(moduleTraitCor)) {
      for(j in 1:ncol(moduleTraitCor)) {
        r_val <- moduleTraitCor[i, j]
        p_val <- moduleTraitPvalue[i, j]
        
        # Format R value (2 decimal places)
        r_text <- sprintf("%.2f", r_val)
        
        # Format P value (scientific notation or standard format)
        if(p_val < 0.001) {
          p_text <- sprintf("(%.1e)", p_val)  # scientific notation, 1 decimal place
        } else if(p_val < 0.01) {
          p_text <- sprintf("(%.3f)", p_val)
        } else {
          p_text <- sprintf("(%.2f)", p_val)
        }
        
        # Combine: R value on new line with P value in parentheses
        textMatrix[i, j] <- paste0(r_text, "\n", p_text)
      }
    }
    
    par(mar = c(8, 10, 4, 4))
    labeledHeatmap(
      Matrix = moduleTraitCor, 
      main = 'Module-trait relationships\n(R: correlation coefficient; P value in parentheses)',
      xLabels = colnames(trait_matched), 
      yLabels = colnames(module_matched), 
      ySymbols = colnames(module_matched),
      colorLabels = FALSE, 
      colors = blueWhiteRed(50), 
      cex.text = 0.5, 
      zlim = c(-1,1),
      textMatrix = textMatrix, 
      setStdMargins = FALSE
    )
    
    dev.off()
    cat("✅ Module-trait correlation heatmap saved: 6_Module_Trait_Heatmap.pdf\n")
    
    # Export correlation results
    write.csv(moduleTraitCor, "moduleTraitCor.csv", row.names = TRUE)
    write.csv(moduleTraitPvalue, "moduleTraitPvalue.csv", row.names = TRUE)
    cat("✅ Correlation results exported:\n")
    cat("   - moduleTraitCor.csv: correlation matrix\n")
    cat("   - moduleTraitPvalue.csv: p-value matrix\n")

    # ============================================================
    # 7. Key module identification: select the most correlated module for each tissue (largest absolute value)
    # ============================================================
    cat("\n=== 7. Key module identification (select the most correlated module for each tissue) ===\n")
    cat("Definition: for each tissue/phenotype, select the module with the largest absolute correlation as the tissue-specific module\n\n")
    
    if(ncol(trait_matched) > 0 && nrow(moduleTraitCor) > 0) {
      
      # Create result summary table
      tissue_specific_modules <- data.frame(
        Trait = character(),
        Top_Module = character(),
        Correlation = numeric(),
        Pvalue = numeric(),
        Gene_Count = integer(),
        stringsAsFactors = FALSE
      )
      
      # Store tissue-specific module information (for Hub gene filtering)
      tissue_module_list <- list()
      
      # Find the most correlated module for each phenotype (tissue)
      for(trait_name in colnames(trait_matched)) {
        cat("\n--- Analyzing phenotype:", trait_name, "---\n")
        
        # Extract correlation coefficients of this phenotype with all modules
        trait_cor <- moduleTraitCor[, trait_name]
        trait_pval <- moduleTraitPvalue[, trait_name]
        
        # Find the module with the largest absolute correlation (tissue-specific module)
        top_idx <- which.max(abs(trait_cor))
        top_module <- names(trait_cor)[top_idx]
        top_cor <- trait_cor[top_idx]
        top_pval <- trait_pval[top_idx]
        
        cat("  ✅ Most correlated module:", top_module, 
            " (r = ", round(top_cor, 3), 
            ", p = ", format(top_pval, scientific = TRUE, digits = 3), ")\n", sep = "")
        
        # Export gene list of this module
        module_color <- sub("^ME", "", top_module)
        gene_module_select <- subset(gene_module, module == module_color)$gene_name
        gene_count <- length(gene_module_select)
        
        if(gene_count > 0) {
          # Create file named after phenotype
          safe_name <- gsub("[^A-Za-z0-9_]", "_", trait_name)
          filename <- paste0("TissueSpecificModule_", safe_name, "_", top_module, ".csv")
          write.csv(data.frame(gene = gene_module_select), 
                    filename, 
                    row.names = FALSE)
          cat("  ✅ Gene list saved:", filename, " (", gene_count, " genes)\n", sep = "")
        } else {
          cat("  ⚠️ Warning: this module has no gene data\n")
        }
        
        # Add to summary table
        tissue_specific_modules <- rbind(tissue_specific_modules, data.frame(
          Trait = trait_name,
          Top_Module = top_module,
          Correlation = round(top_cor, 4),
          Pvalue = format(top_pval, scientific = TRUE, digits = 4),
          Gene_Count = gene_count,
          stringsAsFactors = FALSE
        ))
        
        # Store tissue-specific module information
        tissue_module_list[[trait_name]] <- list(
          module = top_module,
          module_color = module_color,
          correlation = top_cor,
          pvalue = top_pval,
          gene_list = gene_module_select
        )
      }
      
      # ============================================================
      # Export summary of tissue-specific modules
      # ============================================================
      write.csv(tissue_specific_modules, "TissueSpecificModules_summary.csv", row.names = FALSE)
      cat("\n✅ Summary of all tissue-specific modules saved: TissueSpecificModules_summary.csv\n")
      
      cat("\n=== Tissue-specific module summary ===\n")
      print(tissue_specific_modules)
      
      # ============================================================
      # 8. Hub gene filtering (based on tissue-specific modules)
      # ============================================================
      cat("\n=== 8. Hub gene filtering (based on tissue-specific modules) ===\n")
      
      # Set filtering thresholds (adjust as needed)
      MM_threshold <- 0.8      # module membership threshold
      GS_threshold <- 0.2      # gene significance threshold
      p_threshold <- 0.05      # p-value threshold
      
      cat("Filtering criteria:\n")
      cat("  - |MM| >= ", MM_threshold, "\n", sep = "")
      cat("  - |GS| >= ", GS_threshold, "\n", sep = "")
      cat("  - MM_p < ", p_threshold, "\n", sep = "")
      cat("  - GS_p < ", p_threshold, "\n", sep = "")
      
      # Hub gene summary table
      hub_gene_summary <- data.frame(
        Tissue = character(),
        Module = character(),
        Module_Correlation = numeric(),
        Hub_Gene_Count = integer(),
        stringsAsFactors = FALSE
      )
      
      # Store all Hub genes
      all_hub_genes_list <- list()
      
      # Iterate over each tissue and its specific module
      for(trait_name in names(tissue_module_list)) {
        
        cat("\n----------------------------------------\n")
        cat("Analyzing tissue:", trait_name, "\n")
        
        module_info <- tissue_module_list[[trait_name]]
        module_name <- module_info$module
        module_color <- module_info$module_color
        module_cor <- module_info$correlation
        module_pval <- module_info$pvalue
        
        cat("Tissue-specific module:", module_name, "\n")
        cat("Module-tissue correlation:", round(module_cor, 3), 
            " (p = ", format(module_pval, scientific = TRUE, digits = 3), ")\n", sep = "")
        
        # Get gene list of this module
        module_genes <- module_info$gene_list
        
        if(length(module_genes) == 0) {
          cat("⚠️ Module", module_color, "has no genes, skipping\n")
          next
        }
        
        cat("Number of module genes:", length(module_genes), "\n")
        
        # ============================================================
        # Extract expression matrix of this module's genes
        # ============================================================
        module_genes_present <- intersect(module_genes, colnames(gene))
        
        if(length(module_genes_present) == 0) {
          cat("⚠️ Module genes do not exist in expression matrix, skipping\n")
          next
        }
        
        if(length(module_genes_present) < length(module_genes)) {
          cat("⚠️ Warning:", length(module_genes) - length(module_genes_present), 
              "genes do not exist in expression matrix\n")
        }
        
        module_expr <- gene[, module_genes_present, drop = FALSE]
        
        # ============================================================
        # Module Membership (MM) - correlation between gene and module eigengene
        # ============================================================
        
        # Get module eigengene (ME) of this module
        if(!(module_name %in% colnames(module))) {
          cat("⚠️ ", module_name, "not found, skipping\n")
          next
        }
        
        module_ME <- module[, module_name, drop = FALSE]
        
        # Calculate MM
        MM_result <- signedKME(
          module_expr,
          module_ME,
          outputColumnName = "MM"
        )
        
        # Calculate MM p-values
        MM_p <- corPvalueStudent(
          as.matrix(MM_result),
          nrow(module)
        )
        
        # ============================================================
        # Gene Significance (GS) - correlation between gene and corresponding tissue
        # ============================================================
        
        # Extract phenotype data of this tissue
        if(!(trait_name %in% colnames(trait_matched))) {
          cat("⚠️ Phenotype", trait_name, "not found, skipping\n")
          next
        }
        
        tissue_trait <- trait_matched[, trait_name, drop = FALSE]
        
        # Calculate GS
        GS_result <- as.data.frame(
          cor(
            module_expr,
            tissue_trait,
            use = "p"
          )
        )
        
        colnames(GS_result) <- "GS"
        
        # Calculate GS p-values
        GS_p <- corPvalueStudent(
          as.matrix(GS_result),
          nrow(trait_matched)
        )
        
        colnames(GS_p) <- "GS_pvalue"
        
        # ============================================================
        # Merge MM and GS results
        # ============================================================
        
        hub_result <- data.frame(
          Gene = rownames(MM_result),
          MM = MM_result[, 1],
          MM_pvalue = MM_p[, 1],
          GS = GS_result[, 1],
          GS_pvalue = GS_p[, 1],
          stringsAsFactors = FALSE
        )
        
        # ============================================================
        # Hub gene filtering
        # ============================================================
        
        hub_genes_filtered <- subset(
          hub_result,
          abs(MM) >= MM_threshold &
          abs(GS) >= GS_threshold &
          MM_pvalue < p_threshold &
          GS_pvalue < p_threshold
        )
        
        hub_gene_count <- nrow(hub_genes_filtered)
        
        cat("Number of Hub genes:", hub_gene_count, "\n")
        
        # ============================================================
        # Output results
        # ============================================================
        
        if(hub_gene_count > 0) {
          
          # Sort by |GS| descending (larger GS = more important)
          hub_genes_filtered <- hub_genes_filtered[
            order(-abs(hub_genes_filtered$GS)),
          ]
          
          # Save Hub genes for this tissue
          safe_tissue <- gsub("[^A-Za-z0-9_]", "_", trait_name)
          outfile <- paste0(
            "HubGenes_",
            module_color,
            "_specific_for_",
            safe_tissue,
            ".csv"
          )
          
          write.csv(
            hub_genes_filtered,
            outfile,
            row.names = FALSE
          )
          
          cat("✅ Saved:", outfile, "\n")
          
          # Show top 5 Hub genes
          cat("  Top 5 Hub genes:\n")
          print(head(hub_genes_filtered[, c("Gene", "MM", "GS")], 5))
          
          # Add to summary table
          hub_gene_summary <- rbind(
            hub_gene_summary,
            data.frame(
              Tissue = trait_name,
              Module = module_name,
              Module_Correlation = round(module_cor, 4),
              Hub_Gene_Count = hub_gene_count,
              stringsAsFactors = FALSE
            )
          )
          
          # Store in list (for later merging)
          all_hub_genes_list[[paste0(module_color, "_", safe_tissue)]] <- hub_genes_filtered$Gene
          
        } else {
          cat("⚠️ No Hub genes meeting the criteria found\n")
        }
      }
      
      # ============================================================
      # Output Hub gene summary
      # ============================================================
      
      if(nrow(hub_gene_summary) > 0) {
        
        # Sort by Hub gene count descending
        hub_gene_summary <- hub_gene_summary[
          order(-hub_gene_summary$Hub_Gene_Count),
        ]
        
        write.csv(
          hub_gene_summary,
          "Hub_gene_summary.csv",
          row.names = FALSE
        )
        
        cat("\n========================================\n")
        cat("✅ Hub gene filtering completed\n")
        cat("========================================\n")
        
        cat("\n=== Hub gene summary ===\n")
        print(hub_gene_summary)
        
        # ============================================================
        # Optional: output union or intersection of all Hub genes
        # ============================================================
        
        if(length(all_hub_genes_list) > 0) {
          
          # Union of all Hub genes
          all_hub_genes_union <- unique(unlist(all_hub_genes_list))
          write.csv(
            data.frame(Gene = all_hub_genes_union),
            "All_HubGenes_Union.csv",
            row.names = FALSE
          )
          cat("\n✅ Union of all Hub genes saved: All_HubGenes_Union.csv (",
              length(all_hub_genes_union), " genes)\n", sep = "")
          
          # Hub genes appearing in multiple tissues (shared Hub genes)
          gene_counts <- table(unlist(all_hub_genes_list))
          shared_hub_genes <- names(gene_counts[gene_counts >= 2])
          
          if(length(shared_hub_genes) > 0) {
            write.csv(
              data.frame(Gene = shared_hub_genes, 
                         Occurrence = as.integer(gene_counts[shared_hub_genes])),
              "Shared_HubGenes_across_tissues.csv",
              row.names = FALSE
            )
            cat("✅ Shared Hub genes saved: Shared_HubGenes_across_tissues.csv (",
                length(shared_hub_genes), " genes)\n", sep = "")
          }
        }
        
      } else {
        cat("\n⚠️ No Hub genes found\n")
      }
      
    } else {
      cat("⚠️ No available phenotype data or correlation results\n")
    }
    
  }  # ← close else block (common_samples matched successfully)
  
} else {
  cat("⚠️ trait.txt not found, skipping trait association analysis\n")
  cat("Tip: name the phenotype file 'trait.txt' and place it in the current directory\n")
  cat("      Then rerun the script to complete trait association\n")
}  # ← close outermost if(file.exists("trait.txt"))

# ============================================================
# 9. Module eigengene clustering tree and heatmap (inter-module correlation)
# ============================================================
cat("\n=== 9. Module eigengene clustering tree and heatmap ===\n")

if(exists("module") && exists("mergedColors")) {
  
  # Get module eigengenes (MEs)
  # Note: module object is already the module eigengene matrix
  
  cat("Module eigengene dimensions:", dim(module), "\n")
  cat("Number of modules:", ncol(module), "\n")
  
  # ============================================================
  # Plot clustering tree and heatmap of module eigengenes
  # ============================================================
  pdf("7_module_cor_Eigengene-dendrogram.pdf", width = 10, height = 12)
  
  plotEigengeneNetworks(
    module,                       # module eigengene matrix
    setLabels = "",               # leave plot titles empty
    marDendro = c(0, 4, 1, 4),    # margins for dendrogram (bottom, right, top, left)
    marHeatmap = c(3, 4, 1, 2),   # margins for heatmap
    cex.lab = 0.8,                # axis label font size
    xLabelsAngle = 90             # X-axis label rotation angle (90 degrees vertical)
  )
  
  dev.off()
  cat("✅ Module eigengene clustering tree and heatmap saved: 7_module_cor_Eigengene-dendrogram.pdf\n")
    
  # ============================================================
  # Export inter-module correlation matrix
  # ============================================================
  module_cor_matrix <- cor(module, use = 'p')
  write.csv(module_cor_matrix, "module_correlation_matrix.csv", row.names = TRUE)
  cat("✅ Inter-module correlation matrix saved: module_correlation_matrix.csv\n")
  
} else {
  cat("⚠️ Missing module eigengene data, skipping\n")
}

# ============================================================
# 10. Clean up step markers (optional)
# ============================================================
# Uncomment the following code to clean up step markers and force re-run
# cat("\n=== Cleaning up step markers ===\n")
# file.remove(c("step1_done.txt", "step2_done.txt", "step3_done.txt", "step4_done.txt"))

# ============================================================
# 11. Session information
# ============================================================
cat("\n=== Session information ===\n")
cat("R version:", R.version.string, "\n")
cat("WGCNA version:", as.character(packageVersion("WGCNA")), "\n")
cat("Completion time:", date(), "\n")
cat("Working directory:", getwd(), "\n")

cat("\n========================================\n")
cat("     WGCNA Analysis Completed!\n")
cat("========================================\n\n")