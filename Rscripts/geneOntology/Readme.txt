This folder contains scripts for gene ontology enrichment using the R-package topGO. 
There is a functions script and a wrapper script with some example data. 

The example datasets are from TCGA and contain LUAD samples. The file LUAD_PreprocessedData.rda is the expression matrix with all genes,
and LUAD_PreprocessedData_DEgenes.rda is the expression matrix with all differentially expressed genes. The last column in this matrix,
specifies the directionality of the DE gene (up or down). The function for GO-enrichment only needs two vectors of gene names (se comments in functions script) one
contaning all genes (gene univers) and one contaning DE genes (gene of interest) e.g. the counts are not needed. 
Normally, however, the output from the preprocessed data scrips will be the whole count matrix and so the scripts are made for this input.

The first part of the functions script contains code for downloading the human GO-gene associations (which are updated now and again) and making a GO-background
object needed for the enrichement (list object with GO-IDs and associated genes). This part is currently hastagged out and instead a "premade" object is loaded 
directly (HUGO.Rdata). This code can be uncommented to update the HUGO.RData once in a while (this is a heavy computation). 

