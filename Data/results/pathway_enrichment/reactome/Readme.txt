The scripts in this folder may be used to perform the pathway enrichment analysis using the R-package ReactomePA. 
There is a functions script (Reactome_functions.R) and a wrapper script (Reactome_wrapper.R) with an example of how to call the function.

1) Reactome_functions.R contains two functions:
- conversion.background : it converts the genes symbols into 	entrezIDs requested as input by the enrich_pathway function. 
- enrich_pathway : perform the pathway enrichment analysis and a plot is showed as output.

2) Reactome_wrapper.R:
The example datasets are from TCGA and contain LUAD samples. The file LUAD_PreprocessedData.rda is the expression matrix with all genes, and LUAD_PreprocessedData_DEgenes.rda is the expression matrix with all differentially expressed genes. The last column in this matrix, specifies the directionality of the DE gene (up or down).
Both files are located into the folder /data/user/tools_scripts_repository/PATHWAY_2017/example_dataset/.

NB: it could happen that an error occurs after running the enrich_pathway function. Try to increase the "value" argument. I set this value at 0.5 just to get a plot as output, but the 0.05 value is recommended to perform a significant analysis.

