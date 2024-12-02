# Load libraries
library(quarto)
library(fs)

# Create the data
# OK to run just once
# source("data.R") # Creates tse.rds

# Alpha diversity analysis
# Define indices for analysis
indices <- c("shannon", "observed")

# Loop through each index
for (index in indices) {
    # Create temporary qmd with unique name
    temp_qmd <- paste0("alpha_", index, ".qmd")
    file.copy("alpha.qmd", temp_qmd)
    
    # Set parameters for this run
    params_alpha <- list(
        index = index,
        adjust.method = "fdr"
    )
    
    # Render from the temporary file
    quarto_render(
        input = temp_qmd,
        execute_params = params_alpha
    )
    
    # Clean up temporary qmd
    file.remove(temp_qmd)
}

# params_alpha <- list( index = "shannon", 
#                      adjust.method = "fdr")
# quarto_render(input = "alpha.qmd", 
#               execute_params = params_alpha
# )
# 
# # TODO
# quarto::quarto_render("beta.qmd")
# quarto::quarto_render("daa.qmd")

