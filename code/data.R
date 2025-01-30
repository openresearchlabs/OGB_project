library(stringr)
library(readxl)

#load the functions needed
source("funct.R")

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
samdf$gender <- factor(samdf$gender)
samdf$visit <- factor(samdf$visit)
samdf$diet_in <- factor(samdf$diet)
samdf$meal_group <- factor(samdf$meal_group, levels = c("1", "2", "3", "4")) 

# "Meal": intervention Experiment (some hours)
# Half of each diet group received either meal, yielding 4 "diet-meal" groups:
# oats-oats, oats-rice, rice-oats ja rice-rice

samdf <- samdf %>%
  assign_timepoint() %>%
  assign_paired() %>%
  assign_time() %>%
  assign_meal()  %>%
  assign_diet()  # Encode as factors oat/rice instead of 1/2

samdf$timepoint <- factor(samdf$timepoint, levels = c("before", "after"))
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
# Add alpha diversity
tse <- addAlpha(x = tse,assay.type = 'counts',
                index = c('observed', 'shannon'), 
                niter = 100) 

# make primary comparison into tse
# Create a new column for group in colData of the TreeSummarizedExperiment object
# Initialize the group column with NA values directly in tse
tse$group <- NA  

# Assign group labels directly within tse
for (comp in comparisons) {
  # Extract diet and visit information from the comparisons
  group1_info <- extract_diet_visit(comp[1])
  group2_info <- extract_diet_visit(comp[2])
  
  # Assign group labels directly in tse object based on diet and visit
  tse$group[colData(tse)$diet_in == group1_info$diet_in & colData(tse)$visit == group1_info$visit] <- comp[1]
  tse$group[colData(tse)$diet_in == group2_info$diet_in & colData(tse)$visit == group2_info$visit] <- comp[2]
}

# There are some duplicates Id in diet group which prevent the code from 
# performing paired test:
# to remove, and keep first occurrence
# Function to remove duplicates within specific group categories
# Use the function to remove duplicates across both diet groups
tse <- remove_duplicates(tse)
tse <- transformAssay(tse, method = "relabundance")
tse <- agglomerateByRanks(tse)

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
altExp(tse, "phylum") <- agglomerateByVariable(altExp(tse, "phylum"), 
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
altExp(tse, "genus") <- agglomerateByVariable(altExp(tse, "genus"), 
                                   by = "rows", 
                                   f = "genus_sub")

# Agglomerate the data based on specified taxa
altExp(tse, "genus_prevalent") <- agglomerateByRank(tse,  rank="genus")
altExp(tse, "genus_prevalent") <- agglomerateByPrevalence(altExp(tse, "genus_prevalent"), assay.type="relabundance", detection=0.1/100, prevalence=10/100, name="genus_prevalent")

altExp(tse, "phylum_prevalent") <- agglomerateByRank(tse,  rank="phylum")
altExp(tse, "phylum_prevalent") <- agglomerateByPrevalence(altExp(tse, "phylum_prevalent"), assay.type="relabundance", detection=0.1/100, prevalence=10/100, name="phylum_prevalent")

altExp(tse, "family_prevalent") <- agglomerateByRank(tse,  rank="family")
altExp(tse, "family_prevalent") <- agglomerateByPrevalence(altExp(tse, "family_prevalent"), assay.type="relabundance", detection=0.1/100, prevalence=10/100, name="family_prevalent")

altExp(tse, "species_prevalent") <- agglomerateByRank(tse,  rank="species")
altExp(tse, "species_prevalent") <- agglomerateByPrevalence(altExp(tse, "species_prevalent"), assay.type="relabundance", detection=0.1/100, prevalence=10/100, name="species_prevalent")

# Add functional predictions to tse
path_abundance <- read.csv("../data/HUMAnN3/processed/pathabundance_unstratified.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

path_coverage <- read.csv("../data/HUMAnN3/processed/pathcoverage_unstratified.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

genefam_KO <- read.csv("../data/HUMAnN3/final/Renorm_genefamilies_Uniref90_KO_unstratified.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

genefam_metacyc <- read.csv("../data/HUMAnN3/final/Renorm_genefamilies_Uniref90_MetaCyc_unstratified.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)


columns_to_remove <- c("AK1304", "PP2368", "HK2340")
columns_to_keep_pa <- !grepl(paste(columns_to_remove, collapse = "|"), colnames(path_abundance))
columns_to_keep_pc <- !grepl(paste(columns_to_remove, collapse = "|"), colnames(path_coverage))
columns_to_keep_ko <- !grepl(paste(columns_to_remove, collapse = "|"), colnames(genefam_KO))
columns_to_keep_mtc <- !grepl(paste(columns_to_remove, collapse = "|"), colnames(genefam_metacyc))

abundance_matrix_pa <- path_abundance[, columns_to_keep_pa]
abundance_matrix_pc <- path_coverage[, columns_to_keep_pc]
abundance_matrix_ko <- genefam_KO[, columns_to_keep_ko]
abundance_matrix_mtc <- genefam_metacyc[, columns_to_keep_mtc]

colnames(abundance_matrix_pa) <- colnames(tse)
colnames(abundance_matrix_pc) <- colnames(tse)
colnames(abundance_matrix_ko) <- colnames(tse)
colnames(abundance_matrix_mtc) <- colnames(tse)

# Add the filtered matrices to the `AltExp` of the SummarizedExperiment
altExp(tse, "pathabundance") <- SummarizedExperiment(
  assays = list(counts = abundance_matrix_pa),
  rowData = DataFrame(Pathway = rownames(abundance_matrix_pa)),
  colData = colData(tse)  
)

altExp(tse, "pathcoverage") <- SummarizedExperiment(
  assays = list(counts = abundance_matrix_pc),
  rowData = DataFrame(Coverage = rownames(abundance_matrix_pc)),
  colData = colData(tse)  
)

altExp(tse, "KO") <- SummarizedExperiment(
  assays = list(counts = abundance_matrix_ko),
  rowData = DataFrame(Gene_Families = rownames(abundance_matrix_ko)),
  colData = colData(tse)  
)

altExp(tse, "metacyc") <- SummarizedExperiment(
  assays = list(counts = abundance_matrix_mtc),
  rowData = DataFrame(Gene_Families = rownames(abundance_matrix_mtc)),
  colData = colData(tse)  
)

altExp(tse, "pathabundance") <- transformAssay(altExp(tse, "pathabundance"), method = "relabundance")

altExp(tse, "pathcoverage") <- transformAssay(altExp(tse, "pathcoverage"), method = "relabundance")

altExp(tse, "KO") <- transformAssay(altExp(tse, "KO"), method = "relabundance")

altExp(tse, "metacyc") <- transformAssay(altExp(tse, "metacyc"), method = "relabundance")

# Print the group assignments
print(table(tse$group))

# Save TreeSE object for later use
saveRDS(tse, file="../data/tse.Rds")
