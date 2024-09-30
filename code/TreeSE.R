library(mia)
library(TreeSummarizedExperiment)
library(stringr)
library(readxl)
library(dplyr)

# Import metaphlan abundance table as TreeSE object
metaphlan.file <- "../data/metaphlan_db_meta4_combined_reports.txt" # Specify file path
# Import the file as TreeSE
tse <- importMetaPhlAn(metaphlan.file,
                       removeTaxaPrefixes=TRUE
                       ) 

# Simplify the tse colnames and bring them to same format as in sample metadata
colnames(tse) <- gsub("(_.*|\\..*)", "", colnames(tse))

# Import sample metadata
samdf <-  read_excel('../data/metadata.xlsx', sheet = 1, col_names = TRUE)

# Convert to a data frame if necessary
samdf <- as.data.frame(samdf)
rownames(samdf) <- samdf$sample
# Group works better as factor
#samdf$diet <- factor(samdf$diet)
#samdf$visit <- factor(samdf$visit)
samdf$meal_group <- factor(samdf$meal_group)
samdf$gender <- factor(samdf$gender)
# Check that the sample data and assay data match by sample names
if (!all(rownames(samdf)==colnames(tse))) {stop("Check sample ID matching")}

# Add sample metadata to the TreeSE as colData 
colData(tse) <- DataFrame(samdf[colnames(tse),])

#the metaphlan results is essentially relative abundance, so "counts"="relabundance"
#check with colSums(assay(tse, "counts"))
# removing plasmids
tse <- tse[grep("plasmid", rowData(tse)[,"kingdom"],
                                      ignore.case = TRUE, invert = TRUE),]
# Gets a subset of object that includes prevalent taxa, genus level (10% prevalence above 0.1% detection level)
altExp(tse, "PrevalentGenus") <- agglomerateByPrevalence(tse, rank="genus",
                                                             other_label="Other",
                                                             assay.type="counts",
                                                             detection=0.1/100,
                                                             prevalence=10/100)

altExp(tse, "PrevalentSpecies") <- agglomerateByPrevalence(tse, rank="species",
                                                         other_label="Other",
                                                         assay.type="counts",
                                                         detection=0.1/100,
                                                         prevalence=10/100)
tse <- addAlpha(tse, index = c("observed", "shannon"))

#make primary comparison into tse
# Create a new column for group in colData of the TreeSummarizedExperiment object
tse$group <- NA  # Initialize the group column with NA values
sample_data <- as.data.frame(colData(tse))
# Loop through each comparison 
comparisons <- list(
  c("diet_1_visit_1", "diet_1_visit_2"),
  c("diet_2_visit_1", "diet_2_visit_2"),
  c("diet_1_visit_1", "diet_2_visit_1"),
  c("diet_1_visit_2", "diet_2_visit_2")
)

# Loop through each comparison to create plots
for (comp in comparisons) {
  # Parse comparison conditions
  condition_1 <- unlist(strsplit(comp[1], "_"))
  condition_2 <- unlist(strsplit(comp[2], "_"))
  
  # Define group labels
  group1_label <- paste(condition_1, collapse = "_")
  group2_label <- paste(condition_2, collapse = "_")
  
  # Subset sample data for each condition
  group1_samples <- sample_data %>%
    filter(diet == condition_1[2] & visit == condition_1[4]) %>%
    pull(sample)
  
  group2_samples <- sample_data %>%
    filter(diet == condition_2[2] & visit == condition_2[4]) %>%
    pull(sample)
  
  # Print debug information
  cat("Group 1 Samples:", group1_samples, "\n")
  cat("Group 2 Samples:", group2_samples, "\n")
  
  # Check if group1_samples and group2_samples are empty
  if (length(group1_samples) == 0) {
    cat("Warning: No samples found for group1:", group1_label, "\n")
  }
  if (length(group2_samples) == 0) {
    cat("Warning: No samples found for group2:", group2_label, "\n")
  }
  
  # Assign group labels to the corresponding samples in colData
  if (length(group1_samples) > 0) {
    tse$group[colData(tse)$sample %in% group1_samples] <- group1_label
  }
  if (length(group2_samples) > 0) {
    tse$group[colData(tse)$sample %in% group2_samples] <- group2_label
  }
}



# Save TreeSE object for later use
saveRDS(tse, file="../output/tse.Rds")

