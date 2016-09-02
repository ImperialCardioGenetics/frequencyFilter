### Using high-resolution variant frequencies to empower clinical genome interpretation

**Abstract**

Whole exome and genome sequencing have transformed the discovery of genetic variants that cause human Mendelian disease, but discriminating pathogenic from benign variants remains a daunting challenge. Rarity is recognized as a necessary, although not sufficient, criterion for pathogenicity, but frequency cutoffs used in Mendelian analysis are often arbitrary and overly lenient. Recent very large reference datasets, such as the Exome Aggregation Consortium (ExAC), provide an unprecedented opportunity to obtain robust frequency estimates even for very rare variants. Here we present a statistical framework for the frequency-based filtering of candidate disease-causing variants, accounting for disease prevalence, genetic and allelic heterogeneity, inheritance mode, penetrance, and sampling variance in reference datasets. Using the example of cardiomyopathy, we show that our approach reduces by two-thirds the number of candidate variants under consideration in the average exome, and identifies 43 variants previously reported as pathogenic that can now be reclassified. We present precomputed allele frequency cutoffs for all variants in the ExAC dataset.

This repository contains code to reproduce analyses, generate figures, and compile the manuscript.

**Content**

- manuscript -> files to compile the main manuscript and the supplementary material
- src -> source code to compute af_filter values
- data -> output from the src code
- R -> R code to do the analysis and create the figures
- data-raw -> raw data files to recreate the analysis
- figures -> the figures
- ClinVar -> ClinVar submission from variant curation
- packrat -> libraries of packages required by the repository

**Instructions for use with Rstudio**

- Download R studio (https://www.rstudio.com/)
- Install knitr and devtools packages
 `install.packages(c("knitr","devtools"))`
- Download and unpack repository
- Select `File->New project->Existing Directory`
- Browse the 'Project working directory' to match the location of the unpacked repository
- Click `Create Project`
- Wait for the files to be imported and for packrat to install all the required packages (you should get this message: Packrat bootstrap successfully completed. Restarting R and entering packrat mode...)
- Close down Rstudio and restart
- You should now be able to re-create the manuscript and figures by knitting the file: `./manuscript/manuscript.Rmd`
