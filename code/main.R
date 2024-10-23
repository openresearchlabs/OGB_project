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
source("./code/primary2Groups_funct.R")

#load or generate tse
if (file.exists("../data/tse.Rds")) {
  tse <- readRDS("../data/tse.Rds")
} else {
  source("./code/TreeSE.R")
}


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
