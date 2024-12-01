# Load libraries
library(quarto)

# Create the data
# OK to run just once
# source("data.R") # Creates tse.rds

# Alpha diversity analysis
quarto::quarto_render("alpha.qmd") # ../output/alpha.html

# TODO
quarto::quarto_render("beta.qmd")
quarto::quarto_render("daa.qmd")

