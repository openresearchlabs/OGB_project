#the code to run all analysis

# COMPARISONS (PRIMARY)
#•	diet 1, visit 1	vs. 	diet 1, visit 2
#•	diet 2, visit 1	vs. 	diet 2, visit 2
#Did the diet change the microbiota?
#  •	diet 1, visit 1	vs.	diet 2, visit 1
#Was the microbiota equal between the diet groups before treatment?
#  •	diet 1, visit 2	vs.	diet 2, visit 2
#Was there a difference in microbiota between the diet groups after the treatment?

# #RUNNING THE FUNCTIONS
# # List of comparisons
# comparisons <- list(
#   c("diet_1_visit_1", "diet_1_visit_2"),
#   c("diet_2_visit_1", "diet_2_visit_2"),
#   c("diet_1_visit_1", "diet_2_visit_1"),
#   c("diet_1_visit_2", "diet_2_visit_2")
# )

# setwd("code")
library(quarto)
quarto_bin <- quarto::quarto_path()
# quarto::quarto_render("Community_composition.qmd")
# quarto::quarto_render("Alpha_Diversity.qmd")
# quarto::quarto_render("Beta_Diversity.qmd")

# Render the qmd using quarto cli
system(paste(shQuote(quarto_bin), "render Community_composition.qmd --output-dir ../output"))
system(paste(shQuote(quarto_bin), "render Alpha_Diversity.qmd --output-dir ../output"))
system(paste(shQuote(quarto_bin), "render Beta_Diversity.qmd --output-dir ../output"))

