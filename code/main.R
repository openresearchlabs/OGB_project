# Load libraries
library(quarto)

# Create the data
# OK to run just once
# source("data.R") # Creates tse.rds

# Alpha diversity analysis
params_alpha <- list( index = "shannon", 
                     adjust.method = "fdr")
quarto_render(input = "alpha.qmd", 
              execute_params = params_alpha
)
# 
# # TODO
# quarto::quarto_render("beta.qmd")
# quarto::quarto_render("daa.qmd")

