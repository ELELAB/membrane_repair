rm(list=ls())
source("Reactome_functions.R")


#--------------------------------------------------------------------------------------
# example using LUAD dataset : pathway enrichment analysis for up-regulated genes

# NB: recommended p-value cutoff = 0.05
#-----------------------------------------------------------------------------------

#load("/data/user/tools_scripts_repository/PATHWAY_2017/example_dataset/LUAD_PreprocessedData.rda")
#load("/data/user/tools_scripts_repository/PATHWAY_2017/example_dataset/LUAD_PreprocessedData_DEgenes.rda")

#### load data

# time point data used for find the conversion.background gene universe
t3_t0_data <- read.csv("../../results/data_t3_t0.csv", row.names = 1)
t2_t0_data <- read.csv("../../results/data_t2_t0.csv", row.names = 1)
t1_t0_data <- read.csv("../../results/data_t1_t0.csv", row.names = 1)

# DE gene lists
t3_t0_up <- read.csv("../../results/DE_genes/up_data_t3_t0.csv")
t3_t0_down <- read.csv("../../results/DE_genes/down_data_t3_t0.csv")

t2_t0_up <- read.csv("../../results/DE_genes/up_data_t2_t0.csv")
t2_t0_down <- read.csv("../../results/DE_genes/down_data_t2_t0.csv")

t1_t0_up <- read.csv("../../results/DE_genes/up_data_t1_t0.csv")
t1_t0_down <- read.csv("../../results/DE_genes/down_data_t1_t0.csv")





background_t3_t0_data <- convertToPathviewFormat(t3_t0_data)
up_genes_t3_t0 <- as.character(t3_t0_up$geneName)
pathway_up_t3_t0 <- enrich_pathway(background_t3_t0_data,up_genes_t3_t0,"t3_t0_up",TRUE,0.05)

background_t3_t0_data <- convertToPathviewFormat(t3_t0_data)
down_genes_t3_t0 <- as.character(t3_t0_down$geneName)
pathway_down_t3_t0 <- enrich_pathway(background_t3_t0_data,down_genes_t3_t0,"t3_t0_down",TRUE,0.05)


background_t2_t0_data <- convertToPathviewFormat(t2_t0_data)
up_genes_t2_t0 <- as.character(t2_t0_up$geneName)
pathway_up_t2_t0 <- enrich_pathway(background_t2_t0_data,up_genes_t2_t0,"t2_t0_up",TRUE,0.05)

background_t2_t0_data <- convertToPathviewFormat(t2_t0_data)
down_genes_t2_t0 <- as.character(t2_t0_down$geneName)
pathway_down_t2_t0 <- enrich_pathway(background_t2_t0_data,down_genes_t2_t0,"t2_t0_down",TRUE,0.05)


background_t1_t0_data <- convertToPathviewFormat(t1_t0_data)
up_genes_t1_t0 <- as.character(t1_t0_up$geneName)
pathway_up_t1_t0 <- enrich_pathway(background_t1_t0_data,up_genes_t1_t0,"t1_t0_up",TRUE,0.05)

background_t1_t0_data <- convertToPathviewFormat(t1_t0_data)
down_genes_t1_t0 <- as.character(t1_t0_down$geneName)
pathway_down_t1_t0 <- enrich_pathway(background_t1_t0_data,down_genes_t1_t0,"t1_t0_down",TRUE,0.05)
