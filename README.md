# OGB_project

Project collaboration with Prof. Kaisa Linderborg and Enni Mannila.
In this repository you can find the code to [generate TreeSummarizedExperiment object](https://github.com/openresearchlabs/OGB_project/blob/main/code/TreeSE.R) 

The code to run whole analysis can be found in [main.R)](https://github.com/openresearchlabs/OGB_project/blob/main/code/main.R).

To run analysis:

1. Make sure you have a directory called "data" which have the information on taxonomic table, and metadata. We will also save our tse object in that directory. Ideally this will point to the Seafile folder where we put the confidential information.

2. Clone the repository

```git clone https://github.com/openresearchlabs/OGB_project.git```

3. Go to the repository

```cd OGB_project```

2. Run the code
```Rscript ./main.R```


### Notes on the analysis:

For the diversity measure we used the observed species (observed) and shannon alpha diversity index (shannon).
To compare the diversity between 2 group of interest we used wilcoxon.

PCoA plot was also generated based on bray-curtis distance.
And pairwise comparison was performed using PERMANOVA test. 

The species and genus level data was then subsequently analyzed with ANCOM-BC.

All results can be found in the output directory.