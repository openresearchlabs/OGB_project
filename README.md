# Microbiota analysis pipeline for OGB project

This project provides a pipeline for analyzing microbiota data from the OGB (Oat-Gut-Brain) project. It includes alpha diversity, beta diversity, differential abundance analysis, and visualization of results.

The pipeline is implemented as a Quarto website for interactive exploration of results. It processes metagenomic data (after upstream preprocessing) to characterize microbial communities and compare sample groups.

## Repository structure

```
.
├── code/
│   ├── _quarto.yml
│   ├── data.R
│   ├── funct.R
│   ├── main.R
│   └── ...
├── output/               # Rendered HTML reports
├── preprocessing/
│   ├── README.md         # Raw read QC → MetaPhlAn (and optional HUMAnN)
│   ├── raw_reads/        # Puhti / taxprofiler / HUMAnN scripts
│   └── remove_columns_metaphlan_db_meta4_combined_reports.R
└── README.md
```

### Key files

- `code/_quarto.yml` — Quarto website configuration
- `code/main.R` — runs the analysis pipeline
- `code/data.R` — loads and formats data
- `code/funct.R` — analysis functions
- `preprocessing/README.md` — raw read preprocessing (taxprofiler, MetaPhlAn4)
- `preprocessing/remove_columns_metaphlan_db_meta4_combined_reports.R` — cleans MetaPhlAn tables for this repo

## Data flow

```
[ENA raw FASTQ]
    → taxprofiler (fastp, BBduk, human removal) → MetaPhlAn4 profiles
    → optional HUMAnN3 on host-depleted reads
    → remove_columns_metaphlan_*.R
    → code/data.R + Quarto (output/)
```

Upstream steps and scripts: [`preprocessing/README.md`](preprocessing/README.md).

## Usage

### Prerequisites

- R (version 4.0 or higher)
- Quarto (version 1.0 or higher)
- Required R packages (see `code/main.R`)

### Installation

1. Clone this repository:

   ```bash
   git clone https://github.com/openresearchlabs/OGB_project.git
   cd OGB_project
   ```

2. Install required R packages:

   ```r
   install.packages(c("package_name"))  # replace with packages from main.R
   ```

3. Install Quarto: https://quarto.org/docs/get-started/

### Running the analysis

1. Open the project in RStudio or your preferred R environment.

2. Adjust `code/_quarto.yml` if you need different analysis parameters.

3. Run the main script:

   ```r
   source("code/main.R")
   ```

4. View results by opening `output/index.html` in a browser.

### Customizing the analysis

- `code/data.R` — data loading and formatting
- `code/funct.R` — add or change analysis functions
- `code/_quarto.yml` — report structure and content
- `preprocessing/` — upstream metagenomics before this pipeline

## Troubleshooting

- **Memory:** increase the R memory limit if needed:

  ```r
  memory.limit(size=8000)
  ```

- **Missing packages:**

  ```r
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install(c("package1", "package2", "package3"))
  ```

- **Quarto render fails:** check the R console, and confirm Quarto is installed and on your `PATH`.

For other issues, open an issue on the GitHub repository.
