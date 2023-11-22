# BayesDeep
Reconstructing spatial transcriptomics at the single-cell resolution

## Overview
BayesDeep is a novel Bayesian methodology for deeply resolving gene expression for all cells by integrating the molecular profile from spatially resolved transcriptomics (SRT) data and the morphological information extracted from the paired histology image. It builds upon a negative binomial regression model to recover gene expression levels at the single-cell resolution from NGS-based SRT data. 

BayesDeep integrates three distinct modalities from a standard NGS-based SRT experiment: the molecular, image, and geospatial profiles. The molecular profile refers to the spot-resolution gene expression data. The image profile corresponds to the detailed morphological context of the paired histology image in terms of a set of cellular features, which may include cell types, nuclei shape characteristics, and any other relevant explanatory features that can be quantified at large scale. The geospatial profile reveals the spatial relationship between the spots and cells.

The spot-resolution gene expression data, along with the single-cell-resolution morphological features of those cells within spot regions, serve as a reference for recovering the single-cell-resolution gene expressions of all cells, whether within or beyond spot regions. We first modeled the observed read count for a specific gene within a spot using a negative binomial (NB) distribution. Then, the underlying spot-resolution relative gene expression in the NB mean is assumed to be the average of single-cell-resolution relative expression across all cells within the spot. Next, we considered the logarithm of each cell’s relative expression level as a linear combination of covariates. A spike-and-slab prior model is applied for each covariate coefficient. On one hand, this feature selection scheme improves the robustness of our model. On the other hand, the corresponding coefficient matrix B uncovers significant associations between gene expression and cellular characteristics, thereby potentially offering valuable biological insights. With the reconstructed single-cell- resolution spatial molecular profile, we can undertake several pivotal downstream analyses.

![flowchart](flowchart.png)

## Tutorial
For the step-by-step tutorial, please refer to: BayesDeep_tutorial.html

Data used for the analysis can be downloaded from ‘Data_for_BayesDeep’ folder on the Dropbox: https://www.dropbox.com/scl/fo/5wj49vf8561h5nl8mc4pu/h?rlkey=q1xssoiqnenq6ohaxrioflk40&dl=0
