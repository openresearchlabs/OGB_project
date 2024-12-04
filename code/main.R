# Load libraries
library(quarto)
library(fs)

# Create the data
# OK to run just once
# source("data.R") # Creates tse.rds
# Call quarto cli
quarto_bin <- quarto::quarto_path()

# Alpha diversity analysis
# Define indices for analysis
indices <- c("shannon", "observed")
# Loop through each index
for (index in indices) {
    output_dir <- paste0("../output/alpha_", index)
    # Construct the Quarto render command
    render_command <- paste(
        shQuote(quarto_bin), 
        "render alpha.qmd",
        "-P", paste0("index=", index),
        "--output", 
        shQuote(paste0("alpha_", index, ".html")),
        "--output-dir", 
        shQuote(output_dir)
    )
    
    # Execute quarto cli
    system(render_command)
}

# # TODO
# quarto::quarto_render("beta.qmd")
# quarto::quarto_render("daa.qmd")

