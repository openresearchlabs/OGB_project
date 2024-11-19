#the code to run all analysis

#OGB COMPARISONS (PRIMARY)
#•	diet 1, visit 1	vs. 	diet 1, visit 2
#•	diet 2, visit 1	vs. 	diet 2, visit 2
#Did the diet change the microbiota?
#  •	diet 1, visit 1	vs.	diet 2, visit 1
#Was the microbiota equal between the diet groups before treatment?
#  •	diet 1, visit 2	vs.	diet 2, visit 2
#Was there a difference in microbiota between the diet groups after the treatment?

#load the functions needed
source("funct.R")

#load or generate tse
source("TreeSE.R")

#RUNNING THE FUNCTIONS
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

#STEP1: alpha diversity
# Loop through indices and comparisons to create and save plots
final_results <- run_diversity_tests(tse, comparisons, variable, indices, outdir)
# View the final merged results
print(final_results)
for (index in indices) {
  lapply(comparisons, function(comp) create_diversity_plot(tse, comp, index,outdir))
}

#STEP2: PCoA
lapply(comparisons, function(comp) PCoA_plot(tse, comp,variable, outdir))
final_results_permanova <- run_rdba_tests(tse, comparisons, variable, outdir)

#STEP3: ancombc
set.seed(123)
results <- lapply(taxa, function(taxa_level) {
  lapply(comparisons, function(comp) {
    run_ancombc_for_variable(tse, comp, variable, taxa_level)
  })
})

#STEP 4: ancombc mix model 
#not possible: https://github.com/qiime2/q2-composition/issues/133
results <- lapply(taxa, function(taxa_level) {
    run_ancombc_mix(tse,taxa_level)
  })


# STEP 5: Wilcox on taxa
results <- lapply(taxa, function(taxa_level) {
  # Loop over comparisons for the current taxa_level
  lapply(comparisons, function(comp) {
    # Run the Wilcox test
    test_result <- wilcox_test_taxa(tse, comp, variable, taxa_level)
    # Create a file name based on the comparison group
    comp_name <- paste(comp, collapse = "_vs_")
    file_name <- paste0("./output/", "wilcox_test_", taxa_level, "_", comp_name, "_results.csv")
    # Save the result to a CSV file
    write.csv(test_result, file = file_name, row.names = FALSE)
    # Return the file name or test result if needed (optional)
    return(file_name)  # or return(test_result) if you want the result instead
  })
})
