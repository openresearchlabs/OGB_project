#the code to run all analysis

# COMPARISONS (PRIMARY)
#•	diet 1, visit 1	vs. 	diet 1, visit 2
#•	diet 2, visit 1	vs. 	diet 2, visit 2
#Did the diet change the microbiota?
#  •	diet 1, visit 1	vs.	diet 2, visit 1
#Was the microbiota equal between the diet groups before treatment?
#  •	diet 1, visit 2	vs.	diet 2, visit 2
#Was there a difference in microbiota between the diet groups after the treatment?

# Load libraries
library(quarto)

#load the functions needed
source("funct.R")

# Define output directory
outdir   <- "../output/"

# Define the kist of comparisons
comparisons <- list(
  c("diet_1_visit_1", "diet_1_visit_2"),
  c("diet_2_visit_1", "diet_2_visit_2"),
  c("diet_1_visit_1", "diet_2_visit_1"),
  c("diet_1_visit_2", "diet_2_visit_2")
)

# Indices to loop through for alpha diversity plot
indices  <- c("shannon", "observed")
taxa     <- c("genus","species")
variable <- "group"

# Render the qmd in R
quarto::quarto_render("alpha.qmd")
quarto::quarto_render("beta.qmd")
quarto::quarto_render("daa.qmd")

# Option2: render the qmd using quarto cli
# quarto_bin <- quarto::quarto_path() # Set path
#system(paste(shQuote(quarto_bin), "render alpha.qmd --output-dir ../output"))
#system(paste(shQuote(quarto_bin), "render beta.qmd --output-dir ../output"))
#system(paste(shQuote(quarto_bin), "render daa.qmd --output-dir ../output"))

