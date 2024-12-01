library(ANCOMBC)
library(dplyr)
library(DT)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(mia)
library(miaViz)
library(multtest)
library(patchwork)
library(quarto)
library(scater)
library(tidyr)
library(TreeSummarizedExperiment)
library(tidyverse)
library(vegan)


# Define the kist of comparisons
comparisons <- list(
  c("diet_1_visit_1", "diet_1_visit_2"),
  c("diet_2_visit_1", "diet_2_visit_2"),
  c("diet_1_visit_1", "diet_2_visit_1"),
  c("diet_1_visit_2", "diet_2_visit_2")
)



# Define function to extract diet and visit info
extract_diet_visit <- function(comparison) {
  condition <- unlist(strsplit(comparison, "_"))
  list(diet = condition[2], visit = condition[4])
}



run_ancombc_mix <- function(tse,taxa) {
  #extract prevalent
  # Gets a subset of object that includes prevalent taxa, genus level 
  # (10% prevalence above 0.1% detection level)
  altExp(tse, "Prevalent") <- agglomerateByPrevalence(tse, rank=taxa,
                                                      other_label="Other",
                                                      assay.type="relabundance",
                                                      detection=0.1/100,
                                                      prevalence=10/100)
  # Subset the TSE to include only samples for the specified groups in the 
  # comparison                                                   
  #q_col <- paste0("q_", variable, comparison[2]) #q_groupdiet_1_visit_2
  #print(q_col)
  tse_preval <- altExp(tse, "Prevalent")
  # Ensure diet and visit columns are in character or factor form in 
  # colData(tse_preval)
  #colData(tse_preval)$diet <- as.factor(as.character(colData(tse_preval)$diet))
  #colData(tse_preval)$visit <- as.factor(as.character(colData(tse_preval)$visit))
  #colData(tse_preval)$id <- as.factor(as.character(colData(tse_preval)$id))
  #print(colnames(colData(tse_preval)))
  #the interactions diet *visit cannot be achive with tse object, so try to split
  # Extract the abundance data
  abundance_data <- assay(tse_preval, "relabundance")
  # Extract the metadata
  metadata <- as.data.frame(colData(tse_preval))
  # Ensure diet, visit, and id columns exist and are factors
  metadata$diet <- as.factor(as.character(metadata$diet))
  metadata$visit <- as.factor(as.character(metadata$visit))
  metadata$id <- as.factor(as.character(metadata$id))
  #selected assay (prevalent)
  out_taxa <- ancombc2(
    data = abundance_data,
    meta_data=metadata,
    #assay_name = "counts",
    #tax_level = taxa,
    p_adj_method = "fdr",
    prv_cut = 0,
    lib_cut = 0,
    group = "diet",
    fix_formula ="diet * visit",#paste0("diet ", "*" ," visit"),
    rand_formula = "(1 | id)",
    lme_control = lme4::lmerControl(),
    struc_zero = TRUE,
    neg_lb = TRUE,
    alpha = 0.05,
    global = FALSE,
    n_cl = 1,
    verbose = FALSE
  )
  res_taxa <- out_taxa$res
  # # Select columns that contain the variable name
  df_taxa <- res_taxa %>%
    dplyr::select(taxon, contains("diet"))%>%
    dplyr::arrange(names(res_taxa)[ncol(res_taxa) - 4])
  df_taxa_sig <- df_taxa %>% filter(names(df_taxa)[ncol(df_taxa) - 4] < 0.05)
  return(df_taxa_sig)
}

# The following code contain the list of functions that were needed to conduct 
# the 2 groups comparison
# based on primary group of interest

# COMPARISONS (PRIMARY)
#  •	diet 1, visit 1	vs. 	diet 1, visit 2
#  •	diet 2, visit 1	vs. 	diet 2, visit 2
# Did the diet change the microbiota?
#  •	diet 1, visit 1	vs.	diet 2, visit 1
# Was the microbiota equal between the diet groups before treatment?
#  •	diet 1, visit 2	vs.	diet 2, visit 2
# Was there a difference in microbiota between the diet groups after the 
# treatment?


#load tse
#tse_ori <- readRDS("../output/tse.Rds")

# There are some duplicates Id in diet group which prevent the code from 
# performing paired test:
# to remove, and keep first occurrence
# Function to remove duplicates within specific group categories
remove_duplicates <- function(tse) {
  # Extract colData as a data frame for easier manipulation
  df <- as.data.frame(colData(tse))
  # Identify groups (assuming 'group' is the column name for 
  # group categorization)
  unique_groups <- unique(df$group)
  # Initialize a list to store indices of unique samples
  unique_indices <- c()
  # Loop through each group to handle duplicates
  for (group in unique_groups) {
    # Subset data for the current group
    group_data <- df[df$group == group, ]
    # Identify duplicates within the current group
    duplicates <- group_data[duplicated(group_data$id) | duplicated(group_data$id, 
                                                                    fromLast = TRUE), ]
    if (nrow(duplicates) > 0) {
      # Print duplicates found in the current group
      print(paste("Duplicate subjects found in group:", group))
      print(duplicates)
      # Select only the first occurrence of each duplicate within the current 
      # group
      unique_subset <- group_data[!duplicated(group_data$id), ]
      # Alternatively, if you want to select randomly:
      # unique_subset <- group_data %>% slice_sample(n = 1)
      # Store the row indices of unique samples
      unique_indices <- c(unique_indices, rownames(unique_subset))
    } else {
      # If no duplicates, keep all samples from the group
      unique_indices <- c(unique_indices, rownames(group_data))
    }
  }
  
  # Subset the original TSE object to keep only unique samples
  tse_unique <- tse[, unique_indices]
  return(tse_unique)
}

#FUNCTIONS
#STEP1 : richness and alpha diversity
#calculate stats diversity
# Function to perform significance testing and plot
perform_significance_test <- function(tse, comparison, variable, measure) {
  # Subset the TSE to include only samples for the specified groups in the 
  # comparison
  tse_subset <- tse[, colData(tse)$group %in% comparison]
  # Ensure there are no NA values in the variable column
  valid_indices <- complete.cases(colData(tse_subset)[[variable]])
  tse_sub <- tse_subset[, valid_indices]
  # Convert colData to a data frame for easier manipulation
  df <- as.data.frame(colData(tse_sub))
  # Extract diet information for paired/independent comparison decision
  diet_1 <- strsplit(comparison[1], "_")[[1]][2]
  diet_2 <- strsplit(comparison[2], "_")[[1]][2]
  # Determine if test should be paired or independent based on diet
  paired_test <- diet_1 == diet_2 # Check if diet numbers are the same
  # Convert the variable to a factor
  df[[variable]] <- factor(df[[variable]])
  # Extract group data for the Wilcoxon test
  group1 <- df[df[[variable]] == levels(df[[variable]])[1], ]
  group2 <- df[df[[variable]] == levels(df[[variable]])[2], ]
  # print(colnames(group1))
  # Initialize result container for current comparison
  result <- data.frame(variable = variable, measure = measure, 
                       comparison = paste(comparison, collapse = " vs "), 
                       p_value = NA, log2_fold_change = NA, 
                       stringsAsFactors = FALSE)
  # Perform Wilcoxon test only if both groups have enough samples
  if (paired_test) {
    # For paired tests, find common subjects
    common_subjects <- intersect(group1$id, group2$id)
    # print(common_subjects)
    # Filter the groups to include only common subjects
    # Filter the groups to include only common subjects
    group1_paired <- group1[group1$id %in% common_subjects, ]
    group2_paired <- group2[group2$id %in% common_subjects, ]
    # Check if there are enough paired samples
    if (nrow(group1_paired) > 0 & nrow(group2_paired) > 0) {
      # Perform paired Wilcoxon test
      wilcox_test <- wilcox.test(group1_paired[[measure]], 
                                 group2_paired[[measure]], paired = TRUE)
      result$p_value <- wilcox_test$p.value
      # Calculate means and log2 fold change
      mean_group1 <- mean(group1_paired[[measure]], na.rm = TRUE)
      mean_group2 <- mean(group2_paired[[measure]], na.rm = TRUE)
      result$mean_Group1 <- mean_group1
      result$mean_Group2 <- mean_group2
      result$log2_fold_change <- log2(mean_group2 / mean_group1)
      comparison_name <- paste(comparison, collapse = "_vs_")
      result$comparison <- comparison_name
    } else {
      # If no paired samples are found, print a message
      message("No paired samples found for comparison: ", 
              paste(comparison, collapse = " vs "))
    }
  } else {
    # For independent tests, check if both groups have enough samples
    if (nrow(group1) > 0 & nrow(group2) > 0) {
      # Perform independent Wilcoxon test
      wilcox_test <- wilcox.test(group1[[measure]], group2[[measure]], 
                                 paired = FALSE)
      result$p_value <- wilcox_test$p.value
      # Calculate means and log2 fold change
      mean_group1 <- mean(group1[[measure]], na.rm = TRUE)
      mean_group2 <- mean(group2[[measure]], na.rm = TRUE)
      result$mean_Group1 <- mean_group1
      result$mean_Group2 <- mean_group2
      result$log2_fold_change <- log2(mean_group2 / mean_group1)
      comparison_name <- paste(comparison, collapse = "_vs_")
      result$comparison <- comparison_name
    } else {
      # If no samples are found, print a message
      message("No samples found for independent comparison: ", 
              paste(comparison, collapse = " vs "))
    }
  }
  
  return(result)
}



# Function to create and save richness plots for specified comparisons and 
# indices
create_diversity_plot <- function(tse, comparison, index) {
  # Subset the TSE to include only samples for the specified groups in the 
  # comparison
  # tse_subset <- tse[, colData(tse)$group %in% comparison]
  # Create a name for the plot based on the index and comparison
  # Create the richness plot
  diversity_plot <- plotColData(
    tse, 
    y = index, 
    x = "group",
    colour_by = "group",show_boxplot = TRUE, show_violin = FALSE
  ) + 
    geom_signif(comparisons = list(comparison), map_signif_level = FALSE) + 
    theme_bw() + 
    theme(text = element_text(size = 8)) +
    labs(title = paste(index, "Diversity: ", comparison[1], " vs ", 
                       comparison[2])) 
  return(diversity_plot)
}

#STEP2:PCOA
# Perform PCoA
PCoA_plot <- function(tse, comparison, variable) {
    # Subset the TSE to include only samples for the specified groups in the 
    # comparison                                                   
    #tse_subset <- tse[, colData(tse)$group %in% comparison]
  
    tse <- runMDS(
        tse,
        FUN = getDissimilarity,
        method = "bray",
        assay.type = "counts",
        name = "PCoA_BC"
    )
    # Calculate explained variance
    e <- attr(reducedDim(tse, "PCoA_BC"), "eig")
    rel_eig <- e / sum(e[e > 0])

    # Loop over each comparison
    for (comp in comparisons) {
      # Subset the TSE for the specified groups in the comparison
      tse_subset <- tse[, tse$group %in% comp]
      
      # Create the PCoA plot
      p <- plotReducedDim(tse_subset, "PCoA_BC", colour_by = variable)
      
      # Add explained variance for axes
      p1 <- p + labs(
        x = paste("PCoA 1 (", round(100 * rel_eig[[1]], 1), "%)", sep = ""),
        y = paste("PCoA 2 (", round(100 * rel_eig[[2]], 1), "%)", sep = "")
      )
      # Add ellipses for groups
      p_ellipse <- p1 + stat_ellipse(aes(color = tse_subset$group), level = 0.95)
      print(p_ellipse)
    }
}

# Function to perform dbRDA and permanova analysis
perform_dbrda_permanova <- function(tse, variable, comparison) {
    # Subset the TSE to include only samples for the specified groups in the 
    # comparison                                                   
    tse_subset <- tse[, colData(tse)$group %in% comparison]
    # Find indices of rows with complete data in the specified metadata columns
    valid_indices <- complete.cases(colData(tse_subset)[, variable])
    # Subset the tse_transform object to only include rows with complete metadata
    tse_rda <- tse_subset[, valid_indices]
    
    # using mia::getPERMANOVA
    
    perm_res <- getPERMANOVA(tse_rda, assay.type = "relabundance", 
                          formula = X ~ group, 
                          permutations = 999)
    perm_df <- as.data.frame(perm_res)
    return(perm_df)
}

# Function to apply the test over all comparisons and measures
run_rdba_tests <- function(tse, comparisons, variable) {
  # Initialize an empty list to store the results
  all_results <- list()
  # Loop through each comparison
  for (comp in comparisons) {
    # Perform the significance test for the current comparison
    result <- perform_dbrda_permanova(tse, variable, comp)
    # Store result with a meaningful name
    all_results[[paste(comp, collapse = "_vs_")]] <- result  
  }
  return(all_results)
}


#STEP3 : DAA with ANCOM-BC
#ancombc
run_ancombc_for_variable <- function(tse,comparison,variable,taxa) {
  #extract prevalent
  # Gets a subset of object that includes prevalent taxa, genus level 
  # (10% prevalence above 0.1% detection level)
  altExp(tse, "Prevalent") <- agglomerateByPrevalence(tse, rank=taxa,
                                                      other_label="Other",
                                                      assay.type="relabundance",
                                                      detection=0.1/100,
                                                      prevalence=10/100)
  # Subset the TSE to include only samples for the specified groups in the 
  # comparison                                                   
  tse_subset <- tse[, colData(tse)$group %in% comparison]
  # print(paste("Subset groups:", paste(comparison, collapse = ", ")))
  q_col <- paste0("q_", variable, comparison[2]) #q_groupdiet_1_visit_2
  # print(q_col)
  tse_preval <- altExp(tse_subset, "Prevalent")
  # Check if this is a repeated-measures comparison (same diet, different visits)
  same_diet <- strsplit(comparison[1], "_")[[1]][2] == strsplit(comparison[2], "_")[[1]][2]
  rand_formula <- if (same_diet) "(1|id)" else NULL
  #selected assay (prevalent)
  out_taxa <- ancombc2(
    data = tse_preval,
    assay_name = "relabundance",
    tax_level = taxa,
    p_adj_method = "fdr",
    prv_cut = 0,
    lib_cut = 0,
    group = variable,
    fix_formula = paste0(variable),
    rand_formula = rand_formula,
    struc_zero = TRUE,
    neg_lb = TRUE,
    alpha = 0.05,
    global = FALSE,
    lme_control = lme4::lmerControl(),
    n_cl = 1,
    verbose = FALSE
  )
  res_taxa <- out_taxa$res
  # # Select columns that contain the variable name
  df_taxa <- res_taxa %>%
        dplyr::select(taxon, contains(variable))%>%
        dplyr::arrange(!!sym(q_col))
  df_taxa_sig <- df_taxa %>% filter(!!sym(q_col) < 0.05)
  # # Generate appropriate filenames
   comparison_name <- paste(comparison, collapse = "_vs_")
  return(df_taxa_sig)
}


#wilcox test taxa
wilcox_test_taxa <- function(tse, comparison, variable, taxa) {
  altExp(tse, "Prevalent") <- agglomerateByPrevalence(tse, rank=taxa,
                                                      other_label="Other",
                                                      assay.type="counts",
                                                      detection=0.1/100,
                                                      prevalence=10/100)
  # Subset the TSE to include only samples for the specified groups in the 
  # comparison                                                   
  tse_subset <- tse[, colData(tse)$group %in% comparison]
  # Ensure there are no NA values in the variable column
  valid_indices <- complete.cases(colData(tse_subset)[[variable]])
  tse_sub <- tse_subset[, valid_indices]
  tse_preval <- altExp(tse_sub, "Prevalent")
  # Extract the abundance data and sample metadata
  abundance_data <- assay(tse_preval, "counts")
  metadata <- as.data.frame(colData(tse_preval))
  # Ensure the variable is a factor
  metadata[[variable]] <- as.factor(metadata[[variable]])
  # Extract diet information for paired/independent comparison decision
  diet_1 <- strsplit(comparison[1], "_")[[1]][2]
  diet_2 <- strsplit(comparison[2], "_")[[1]][2]
  # Determine if test should be paired or independent based on diet
  paired_test <- diet_1 == diet_2 # Check if diet numbers are the same
  # Extract group data for the Wilcoxon test
  group1 <- metadata[metadata[[variable]] == levels(metadata[[variable]])[1], ]
  group2 <- metadata[metadata[[variable]] == levels(metadata[[variable]])[2], ]
  # Initialize a list to store results
  test_results <- list()
  for (taxon in rownames(abundance_data)) {
  # Initialize result container for current comparison
  taxon_result <- data.frame(taxon = taxon, variable = variable, stringsAsFactors = FALSE)
  # Perform Wilcoxon test only if both groups have enough samples
  if (paired_test) {
    # For paired tests, find common subjects
    common_subjects <- intersect(group1$id, group2$id)
    print(common_subjects)
    # Filter the groups to include only common subjects
    group1_paired <- group1[group1$id %in% common_subjects, ]
    group2_paired <- group2[group2$id %in% common_subjects, ]
    # Check if there are enough paired samples
    if (nrow(group1_paired) > 0 & nrow(group2_paired) > 0) {
      # Perform paired Wilcoxon test
      wilcox_test <- wilcox.test(abundance_data[taxon, group1_paired$sample],
                                   abundance_data[taxon, group2_paired$sample],
                                   paired = TRUE)
      taxon_result$p_value <- wilcox_test$p.value
      # Calculate means and log2 fold change
      mean_group1 <- mean(abundance_data[taxon, group1_paired$sample], na.rm = TRUE)
      mean_group2 <- mean(abundance_data[taxon, group2_paired$sample], na.rm = TRUE)
      taxon_result$mean_Group1 <- mean_group1
      taxon_result$mean_Group2 <- mean_group2
      taxon_result$log2_fold_change <- log2(mean_group2 / mean_group1)
      comparison_name <- paste(comparison, collapse = "_vs_")
      taxon_result$comparison <- comparison_name
    } else {
      # If no paired samples are found, print a message
      message("No paired samples found for comparison: ", paste(comparison, collapse = " vs "))
    }
  } else {
    # For independent tests, check if both groups have enough samples
    if (nrow(group1) > 0 & nrow(group2) > 0) {
      # Perform independent Wilcoxon test
      wilcox_test <- wilcox.test(abundance_data[taxon, group1$sample],
                                   abundance_data[taxon, group2$sample],
                                   paired = FALSE)
      taxon_result$p_value <- wilcox_test$p.value
      # Calculate means and log2 fold change
      mean_group1 <- mean(abundance_data[taxon, group1$sample], na.rm = TRUE)
      mean_group2 <- mean(abundance_data[taxon, group2$sample], na.rm = TRUE)
      taxon_result$mean_Group1 <- mean_group1
      taxon_result$mean_Group2 <- mean_group2
      taxon_result$log2_fold_change <- log2(mean_group2 / mean_group1)
      comparison_name <- paste(comparison, collapse = "_vs_")
      taxon_result$comparison <- comparison_name
    } else {
      # If no samples are found, print a message
      message("No samples found for independent comparison: ", paste(comparison, collapse = " vs "))
    }
  }
    # Store the result for the current taxon
    test_results[[taxon]] <- taxon_result
  }
    # Combine results into a single data frame
  test_results_df <- bind_rows(test_results)
  # Adjust p-values for multiple testing (FDR)
  test_results_df$p_adjusted <- p.adjust(test_results_df$p_value, method = "fdr")
  # Sort the results by adjusted p-value
  test_results_df <- test_results_df %>%
    arrange(p_adjusted)
  return(test_results_df)
}
