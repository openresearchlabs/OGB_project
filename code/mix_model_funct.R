library(scater)
library(ggplot2)
library(ggsignif)
library(tidyr)
library(multtest)
library(ANCOMBC)
library(vegan)

run_ancombc_mix <- function(tse,taxa) {
  #extract prevalent
  # Gets a subset of object that includes prevalent taxa, genus level (10% prevalence above 0.1% detection level)
  altExp(tse, "Prevalent") <- agglomerateByPrevalence(tse, rank=taxa,
                                                      other_label="Other",
                                                      assay.type="counts",
                                                      detection=0.1/100,
                                                      prevalence=10/100)
  # Subset the TSE to include only samples for the specified groups in the comparison                                                   
  #q_col <- paste0("q_", variable, comparison[2]) #q_groupdiet_1_visit_2
  #print(q_col)
  tse_preval <- altExp(tse, "Prevalent")
  # Ensure diet and visit columns are in character or factor form in colData(tse_preval)
  #colData(tse_preval)$diet <- as.factor(as.character(colData(tse_preval)$diet))
  #colData(tse_preval)$visit <- as.factor(as.character(colData(tse_preval)$visit))
  #colData(tse_preval)$id <- as.factor(as.character(colData(tse_preval)$id))
  #print(colnames(colData(tse_preval)))
  #the interactions diet *visit cannot be achive with tse object, so try to split
  # Extract the abundance data
  abundance_data <- assay(tse_preval, "counts")
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
    verbose = TRUE
  )
  res_taxa <- out_taxa$res
  # # Select columns that contain the variable name
  df_taxa <- res_taxa %>%
    dplyr::select(taxon, contains(diet))%>%
    dplyr::arrange(!!as.name(names(res_taxa)[ncol(res_taxa) - 1]))
  df_taxa_sig <- df_taxa %>% filter(!!as.name(names(df_taxa)[ncol(df_taxa) - 1]) < 0.05)
  # 
  # # Save results
  write.csv(df_taxa, file = paste0(outdir, "mix_ancombc_", taxa, "_results_", ".csv"), row.names = FALSE)
  write.csv(df_taxa_sig, file = paste0(outdir, "mix_significant_ancombc_", taxa, "_results_", ".csv"), row.names = FALSE)
  saveRDS(res_taxa, file = paste0(outdir, "mix_ancombc_", taxa, "_results_", ".rds"))
}