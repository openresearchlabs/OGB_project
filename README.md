# OGB_project

Project collaboration with Prof. Kaisa Linderborg and Enni Mannila.
In this repository you can find the code to [generate TreeSummarizedExperiment object](https://github.com/openresearchlabs/OGB_project/blob/main/code/TreeSE.R) 

The code to run whole analysis can be found in [main.R)](https://github.com/openresearchlabs/OGB_project/blob/main/code/main.R).

To run analysis:

1. Make sure you have a directory called "data" which have the information on taxonomic table, and metadata. We will also save our tse object in that directory. Ideally this will point to the Seafile folder where we put the confidential information.

2. Clone the repository

```bash
git clone https://github.com/openresearchlabs/OGB_project.git
```

3. Go to the repository

```bash
cd OGB_project/
```
4. Navigate to the directory

```bash
cd code/
```

5. Run the code

```bash
Rscript "main.R"
```
6. Run the Quarto document for alpha diversity in RStudio

```bash
setwd("code")
library(quarto)
quarto::quarto_render("Community_composition_OGB")
quarto::quarto_render("Alpha_Diversity_OGB")
quarto::quarto_render("Beta_Diversity_OGB")
```
