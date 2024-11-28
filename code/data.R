library(stringr)
library(readxl)

# Indices to loop through for alpha diversity plot
indices <- c("shannon", "observed")
taxa <- c("genus","species")
outdir ="./output/"
variable <- "group"

source("funct.R")

# # STEP 5: Wilcox on taxa
# results <- lapply(taxa, function(taxa_level) {
#   # Loop over comparisons for the current taxa_level
#   lapply(comparisons, function(comp) {
#     # Run the Wilcox test
#     test_result <- wilcox_test_taxa(tse, comp, variable, taxa_level)
#     # Create a file name based on the comparison group
#     comp_name <- paste(comp, collapse = "_vs_")
#     file_name <- paste0("./output/", "wilcox_test_", taxa_level, "_", comp_name, "_results.csv")
#     # Save the result to a CSV file
#     write.csv(test_result, file = file_name, row.names = FALSE)
#     # Return the file name or test result if needed (optional)
#     return(file_name)  # or return(test_result) if you want the result instead
#   })
# })


# Import metaphlan abundance table as TreeSE object
metaphlan.file <- "../data/modified_metaphlan_db_meta4_combined_reports.txt" 
# Specify file path
# Import the file as TreeSE
tse <- importMetaPhlAn(metaphlan.file,
                       removeTaxaPrefixes=TRUE
) 
# Simplify the tse colnames and bring them to same format as in sample metadata
colnames(tse) <- gsub("(_.*|\\..*)", "", colnames(tse))

# Import sample metadata
samdf <-  read_excel("../data/metadata.xlsx", sheet = 1, col_names = TRUE)

# Convert to a data frame if necessary
samdf <- as.data.frame(samdf)
rownames(samdf) <- samdf$sample
# Group works better as factor
samdf$diet <- factor(samdf$diet)
samdf$visit <- factor(samdf$visit)
samdf$meal_group <- factor(samdf$meal_group)
samdf$gender <- factor(samdf$gender)
# Check that the sample data and assay data match by sample names
if (!all(rownames(samdf)==colnames(tse))) {stop("Check sample ID matching")}

# Add sample metadata to the TreeSE as colData 
colData(tse) <- DataFrame(samdf[colnames(tse),])

#the metaphlan results is essentially relative abundance, so "counts"="relabundance"
#check with colSums(assay(tse, "counts"))
colSums(assay(tse, "counts"))
# removing plasmids
tse <- tse[grep("plasmid", rowData(tse)[,"kingdom"],
                ignore.case = TRUE, invert = TRUE),]
# Gets a subset of object that includes prevalent taxa, genus level (10% prevalence above 0.1% detection level)
altExp(tse, "PrevalentGenus") <- agglomerateByPrevalence(tse, rank="genus",
                                                         other_label="Other",
                                                         assay.type="counts",
                                                         detection=0.1/100,
                                                         prevalence=10/100,)

altExp(tse, "PrevalentSpecies") <- agglomerateByPrevalence(tse, rank="species",
                                                           other_label="Other",
                                                           assay.type="counts",
                                                           detection=0.1/100,
                                                           prevalence=10/100)

tse <- addAlpha(tse, index = c("observed", "shannon"))

#make primary comparison into tse
# Create a new column for group in colData of the TreeSummarizedExperiment object
# Initialize the group column with NA values directly in tse
tse$group <- NA  

# Define the comparison groups
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

# Assign group labels directly within tse
for (comp in comparisons) {
  # Extract diet and visit information from the comparisons
  group1_info <- extract_diet_visit(comp[1])
  group2_info <- extract_diet_visit(comp[2])
  
  # Assign group labels directly in tse object based on diet and visit
  tse$group[colData(tse)$diet == group1_info$diet & colData(tse)$visit == group1_info$visit] <- comp[1]
  tse$group[colData(tse)$diet == group2_info$diet & colData(tse)$visit == group2_info$visit] <- comp[2]
}

# There are some duplicates Id in diet group which prevent the code from 
# performing paired test:
# to remove, and keep first occurrence
# Function to remove duplicates within specific group categories
# Use the function to remove duplicates across both diet groups
tse <- remove_duplicates(tse)
tse <- transformAssay(tse, method = "relabundance")
tse <- agglomerateByRanks(tse)
tse <- addAlpha(x = tse,assay.type = 'counts',
                index = c('observed', 'shannon', 'simpson'), 
                niter = 100) 
# Changes old levels with new levels
tse$group <- factor(tse$group)

# Getting top taxa on a Phylum level
top_phyla <- getTop(altExp(tse, "phylum"), top = 4, assay.type = "relabundance")

# Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
phylum_renamed <- lapply(rowData(altExp(tse, "phylum"))$phylum, function(x) {
  if (x %in% top_phyla) { x } else { "Other" }
})
rowData(altExp(tse, "phylum"))$phylum_sub <- as.character(phylum_renamed)

# Agglomerate the data based on specified phyla
tse_phylum <- agglomerateByVariable(altExp(tse, "phylum"), 
                                    by = "rows", 
                                    f = "phylum_sub")

# Getting top taxa on a Genus level
top_genera <- getTop(altExp(tse, "genus"), top = 10, assay.type = "relabundance")

# Renaming the "Genus" rank to keep only top taxa and the rest to "Other"
genus_renamed <- lapply(rowData(altExp(tse, "genus"))$genus, function(x) {
  if (x %in% top_genera) { x } else { "Other" }
})
rowData(altExp(tse, "genus"))$genus_sub <- as.character(genus_renamed)

# Agglomerate the data based on specified taxa
tse_genus <- agglomerateByVariable(altExp(tse, "genus"), 
                                   by = "rows", 
                                   f = "genus_sub")

# Print the group assignments
print(table(tse$group))

# Save TreeSE object for later use
# saveRDS(tse, file="../data/tse.Rds")
