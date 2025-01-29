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
lapply(indices, function(index) {
    output_dir <- path(path_abs("../output/alpha"), paste0("alpha_", index))

    # HTML render
    system(paste(
        shQuote(quarto_bin),
        "render alpha/alpha.qmd",
        "-P", paste0("index=", index),
        "--to html",
        "--output", shQuote(paste0("alpha_", index, ".html")),
        "--output-dir", shQuote(output_dir)
    ))

    # PDF render
    system(paste(
        shQuote(quarto_bin),
        "render alpha/alpha.qmd",
        "-P", paste0("index=", index),
        "--to pdf",
        "--output", shQuote(paste0("alpha_", index, ".pdf")),
        "--output-dir", shQuote(output_dir)
    ))
})
# # Beta diversity analysis
system(paste(shQuote(quarto_bin), "render beta/beta.qmd --output-dir ../output/beta"))
# # B/F ratio
system(paste(shQuote(quarto_bin), "render ratio/ratio.qmd --output-dir ../output/ratio"))

# DAA
system(paste(shQuote(quarto_bin), "render daa/daa_taxa.qmd --output-dir ../output"))
system(paste(shQuote(quarto_bin), "render daa/daa_pa.qmd --output-dir ../output"))
system(paste(shQuote(quarto_bin), "render daa/daa_pc.qmd --output-dir ../output"))
system(paste(shQuote(quarto_bin), "render daa/daa_ko.qmd --output-dir ../output"))
system(paste(shQuote(quarto_bin), "render daa/daa_mtc.qmd --output-dir ../output"))

