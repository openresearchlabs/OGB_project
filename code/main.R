# Load libraries
library(quarto)
library(fs)
library(rlang)

# Create the data
# OK to run just once
# source("data.R") # Creates tse.rds

# Alpha diversity analysis
indices <- c("shannon", "observed")

# Render alpha diversity with different parameters
lapply(indices, function(index) {
    orig_dir <- dirname("alpha/alpha.qmd")
    temp_qmd <- file.path(orig_dir, paste0("alpha_", index, ".qmd"))
    file.copy("alpha/alpha.qmd", temp_qmd)
    
    quarto::quarto_render(
        input = temp_qmd,
        execute_params = list(index = index)
    )
    file.remove(temp_qmd)
})

# Finally, render the entire website for qmds that do not directly take params from mainR
quarto::quarto_render()
