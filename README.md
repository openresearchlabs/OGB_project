# OGB_project

Project collaboration with Prof. Kaisa Linderborg and Enni Mannila.
In this repository you can find the code to [generate TreeSummarizedExperiment object](https://github.com/openresearchlabs/OGB_project/blob/main/code/TreeSE.R)

Raw data was preprocessed for all the samples with correct annotations. The output files are available in Seafile data directory.

In the provided data and metadata files, some adjustments were made as follows:
Sample exclusions:
The following {sample(id)} were excluded from the analysis as they were either dropouts or unable to tolerate the study diet, and therefore were not required for the analyses:
AK1304 (304AK), PP2368 (368PP), and HK2340 (340HK).

Visit category annotations:
Annotations for the visit categories of the following {sample(id)} were corrected:
JK1312 (312JK), HP1337 (337HP), AM2310 (310AM), and AM1314 (314AM).

Reference for Visit Annotations:
The “1” in the ID code (e.g., XX1300) refers to the baseline visit,
The “2” in the ID code (e.g., XX2300) refers to the 6th-week visit.

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
