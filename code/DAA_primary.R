#The code to perform
#Differential analysis 

library(mia)
library(dplyr)
library(scater)
library(ggplot2)
library(ggsignif)
library(tidyr)
library(multtest)
library(ANCOMBC)

#OGB COMPARISONS (PRIMARY)
#•	diet 1, visit 1	vs. 	diet 1, visit 2
#•	diet 2, visit 1	vs. 	diet 2, visit 2
#Did the diet change the microbiota?
#  •	diet 1, visit 1	vs.	diet 2, visit 1
#Was the microbiota equal between the diet groups before treatment?
#  •	diet 1, visit 2	vs.	diet 2, visit 2
#Was there a difference in microbiota between the diet groups after the treatment?

#load tse
tse <- readRDS("../output/tse.Rds")

#alpha diversity
#calculate stats diversity
# Function to perform significance testing and plot
perform_significance_test <- function(tse, comparison, variable, measure) {
  # Subset the TSE to include only samples for the specified groups in the comparison
  tse_subset <- tse[, colData(tse)$group %in% comparison]
  # Ensure there are no NA values in the variable column
  valid_indices <- complete.cases(colData(tse_subset)[[variable]])
  tse_sub <- tse_subset[, valid_indices]
  # Convert colData to a data frame for easier manipulation
  df <- as.data.frame(colData(tse_sub))
  # Convert the variable to a factor
  df[[variable]] <- factor(df[[variable]])
  # Extract group data for the Wilcoxon test
  group1 <- df[df[[variable]] == levels(df[[variable]])[1], measure]
  group2 <- df[df[[variable]] == levels(df[[variable]])[2], measure]
  # Initialize result container for current comparison
  result <- data.frame(variable = variable, measure = measure, comparison = paste(comparison, collapse = " vs "), p_value = NA, log2_fold_change = NA, stringsAsFactors = FALSE)
  # Perform Wilcoxon test only if both groups have enough samples
  if (length(group1) > 1 & length(group2) > 1) {
    wilcox_test <- wilcox.test(group1, group2)
    result$p_value <- wilcox_test$p.value
    # Calculate means and log2 fold change
    mean_group1 <- mean(group1, na.rm = TRUE)
    mean_group2 <- mean(group2, na.rm = TRUE)
    result$mean_Group1 <- mean_group1
    result$mean_Group2 <- mean_group2
    result$log2_fold_change <- log2(mean_group2 / mean_group1)

    comparison_name <- paste(comparison, collapse = "_vs_")
    result$comparison <- comparison_name
  }
  return(result)
}

# Function to apply the test over all comparisons and measures
run_diversity_tests <- function(tse, comparisons, variable, indices, outdir) {
  # Initialize an empty list to store the results
  all_results <- list()
  # Loop through each diversity index (e.g., shannon, observed)
  for (measure in indices) {
    # Loop through each comparison
    comparison_results <- lapply(comparisons, function(comp) {
      # Perform the significance test for the current comparison
      result <- perform_significance_test(tse, comp, variable, measure)
      return(result)
    })
    # Combine results for this measure across all comparisons
    comparison_results_df <- do.call(rbind, comparison_results)
    # Add the current index's result to the overall list
    all_results[[measure]] <- comparison_results_df
  }
  # Combine all the results into one data frame
  final_results <- do.call(rbind, all_results)
  # Adjust p-values for multiple testing
  final_results$p_adjusted <- p.adjust(final_results$p_value, method = "holm")
  # Save the final merged table to a CSV file
  write.csv(final_results, file = paste0(outdir, "merged_diversity_results.csv"), row.names = FALSE)
  return(final_results)
}


# Function to create and save richness plots for specified comparisons and indices
create_diversity_plot <- function(tse, comparison, index,outdir) {
  # Subset the TSE to include only samples for the specified groups in the comparison
  tse_subset <- tse[, colData(tse)$group %in% comparison]
  
  # Create a name for the plot based on the index and comparison
  plot_name <- paste0(outdir,index, "_", comparison[1], "_vs_", comparison[2], ".pdf")
  
  # Create the richness plot
  diversity_plot <- plotColData(
    tse_subset, 
    y = index, 
    x = "group",
    colour_by = "group",show_boxplot = TRUE, show_violin = FALSE
  ) + 
    geom_signif(comparisons = list(comparison), map_signif_level = FALSE) + 
    theme_bw() + 
    theme(text = element_text(size = 8)) +
    labs(title = paste(index, "Diversity: ", comparison[1], " vs ", comparison[2]))  
  
  # Save the plot as PDF
  ggsave(filename = plot_name, plot = diversity_plot, width = 8, height = 6, units = "in")
}


#ancombc
run_ancombc_for_variable <- function(tse,comparison,variable,taxa) {
#extract prevalent
# Gets a subset of object that includes prevalent taxa, genus level (10% prevalence above 0.1% detection level)
  altExp(tse, "Prevalent") <- agglomerateByPrevalence(tse, rank=taxa,
                                                      other_label="Other",
                                                      assay.type="counts",
                                                      detection=0.1/100,
                                                      prevalence=10/100)
  # Subset the TSE to include only samples for the specified groups in the comparison                                                   
  tse_subset <- tse[, colData(tse)$group %in% comparison]
  print(paste("Subset groups:", paste(comparison, collapse = ", ")))
  q_col <- paste0("q_", variable, comparison[2]) #q_groupdiet_1_visit_2
  print(q_col)
  tse_preval <- altExp(tse_subset, "Prevalent")
  #selected assay (prevalent)
  out_taxa <- ancombc2(
    data = tse_preval,
    assay_name = "counts",
    tax_level = taxa,
    p_adj_method = "fdr",
    prv_cut = 0,
    lib_cut = 0,
    group = variable,
    fix_formula = variable,
    struc_zero = TRUE,
    neg_lb = TRUE,
    alpha = 0.05,
    global = FALSE,
    n_cl = 1,
    verbose = TRUE
  )

  res_taxa <- out_taxa$res
  # # Select columns that contain the variable name
  df_taxa <- res_taxa %>%
        dplyr::select(taxon, contains(variable))%>%
        dplyr::arrange(!!sym(q_col))
  df_taxa_sig <- df_taxa %>% filter(!!sym(q_col) < 0.05)
  # # Generate appropriate filenames
   comparison_name <- paste(comparison, collapse = "_vs_")
  # 
  # # Save results
   write.csv(df_taxa, file = paste0(outdir, "ancombc_", taxa, "_results_", comparison_name, ".csv"), row.names = FALSE)
   write.csv(df_taxa_sig, file = paste0(outdir, "significant_ancombc_", taxa, "_results_", comparison_name, ".csv"), row.names = FALSE)
   saveRDS(res_taxa, file = paste0(outdir, "ancombc_", taxa, "_results_", comparison_name, ".rds"))
}


#PCOA
# Perform PCoA
PCoA_plot <- function(tse, comparison, variable, outdir) {
  # Subset the TSE to include only samples for the specified groups in the comparison                                                   
tse_subset <- tse[, colData(tse)$group %in% comparison]
tse_subset <- runMDS(
  tse_subset,
  FUN = getDissimilarity,
  method = "bray",
  assay.type = "counts",
  name = "MDS_bray"
)

# Create ggplot object
p <- plotReducedDim(tse_subset, "MDS_bray", colour_by = variable)

# Calculate explained variance
e <- attr(reducedDim(tse_subset, "MDS_bray"), "eig")
rel_eig <- e / sum(e[e > 0])

# Add explained variance for each axis
p1 <- p + labs(
  x = paste("PCoA 1 (", round(100 * rel_eig[[1]], 1), "%", ")", sep = ""),
  y = paste("PCoA 2 (", round(100 * rel_eig[[2]], 1), "%", ")", sep = "")
)
p1
p_ellipse <- p1 + stat_ellipse(aes(color = colour_by), level = 0.95)
#Save the plot as PDF
plot_name <- paste0(outdir,"PCoA", "_", comparison[1], "_vs_", comparison[2], ".pdf")
ggsave(filename = plot_name, plot = p_ellipse, width = 8, height = 6, units = "in")
}
# List of comparisons
comparisons <- list(
  c("diet_1_visit_1", "diet_1_visit_2"),
  c("diet_2_visit_1", "diet_2_visit_2"),
  c("diet_1_visit_1", "diet_2_visit_1"),
  c("diet_1_visit_2", "diet_2_visit_2")
)
# Indices to loop through for alpha diversity plot
indices <- c("shannon", "observed")
taxa <- c("genus","species")
outdir ="./output/"
variable <- "group"

#alpha diversity
# Loop through indices and comparisons to create and save plots
final_results <- run_diversity_tests(tse, comparisons, variable, indices, outdir)
# View the final merged results
print(final_results)

for (index in indices) {
  lapply(comparisons, function(comp) create_diversity_plot(tse, comp, index,outdir))
}

#PCoA
lapply(comparisons, function(comp) PCoA_plot(tse, comp,variable, outdir))


#ancombc
set.seed(123)
results <- lapply(taxa, function(taxa_level) {
  lapply(comparisons, function(comp) {
    run_ancombc_for_variable(tse, comp, variable, taxa_level)
  })
})

