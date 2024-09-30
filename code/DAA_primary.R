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
for (index in indices) {
  lapply(comparisons, function(comp) create_diversity_plot(tse, comp, index,outdir))
}

#ancombc
set.seed(123)
results <- lapply(taxa, function(taxa_level) {
  lapply(comparisons, function(comp) {
    run_ancombc_for_variable(tse, comp, variable, taxa_level)
  })
})

