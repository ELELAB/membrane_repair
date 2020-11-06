library(GOplot)
source("geneOntology_functions.R")

load("HUGO.RData")

#--------------------------------------------------------------------------------------
#Function that  makes genes dataframe to input in the circle_dat function

# fileName: name of file that includes logFC. It is usually a DEA output file
# genes_interest: a character vector that includes all the genes that belong to GO terms
# genes (output): A data frame with columns for 'ID' (genes symbol), 'logFC'

#--------------------------------------------------------------------------------------
make_genes_dataframe <- function(fileName,genes_interest){
  genes <- read.csv(fileName)
  ID <- rownames(genes)
  logFC <- genes$logFC
  genes <- data.frame(ID,logFC)
  genes <- subset(genes,genes$ID %in% genes_interest)
  rownames(genes) <- 1:nrow(genes)
  return(genes)
}


make_genes_column <- function(genes){ 
  c <- genes[1]
  if(length(genes)>=2){
    for(i in 2:length(genes)){
      c <- paste(c,genes[i],sep=",")
      c <- c
    }
  }
  return(c)
}

make_vector <- function(GO.up, my.up.genes, HUGO){
  vector <- c()
  for(i in 1:nrow(GO.up)){
    genes <- intersect(my.up.genes,HUGO[[GO.up$GO.ID[i]]])
    if(length(genes)==0)
      c <- 0
    else{
      c <- make_genes_column(genes)
  }
  vector <- c(vector,c)
  }
  return(vector)
}

make_genes_interest <- function(GO.up, my.up.genes, HUGO){
  genes_interest <- c()
  GO <- as.character(GO.up$GO.ID)
  for(i in 1:nrow(GO.up)){
    genes <- intersect(my.up.genes,HUGO[[GO[i]]])
    genes_interest <- c(genes_interest,genes)
  }
  genes_interest <- genes_interest[!duplicated(genes_interest)]
  return(genes_interest)
}

#------------------------------------------------------------------------------------
# make terms dataframe to input in the circle_dat function
#-------------------------------------------------------------------------------------
make_terms_dataframe <- function(GO.up,vector, category){
  category <- category
  ID <- GO.up$GO.ID
  term <- GO.up$Term
  adj_pval <- GO.up$FDR_fisher
  genes <- vector
  terms <- data.frame(category,ID,term,adj_pval,genes)
  # remove the terms that not include any genes
  terms <- subset(terms,!terms$genes==0)
  return(terms)
}

#-------------------------------------------------------------------------------
# Function that makes circlize plot
# terms: dataframe with 4 columns (category,ID,term,adj_pvalue,genes)
# genes: A data frame with columns for 'ID', 'logFC'
# plot: plot name
#-----------------------------------------------------------------------------
circle_plot <- function(terms,genes,plot){
 
  circ <- circle_dat(terms, genes)
  #process: A character vector of selected processes
  process <- terms$term
  chord <- chord_dat(data = circ, genes, process)
  
  #png(filename=plot,height = 1000, width = 1500)
  #GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 4,lfc.min = min(circ$logFC), lfc.max = max(circ$logFC))
  #dev.off()
  return(chord)
}




# GO terms made by TOPGO (BP)
t1_t0_GO_down <- read.csv("../../Data/results/GO/BP/t1_t0_GO_down.csv", row.names = 1, as.is = TRUE)
t1_t0_GO_up <- read.csv("../../Data/results/GO/BP/t1_t0_GO_up.csv", row.names = 1, as.is = TRUE)

t2_t0_GO_down <- read.csv("../../Data/results/GO/BP/t2_t0_GO_down.csv", row.names = 1, as.is = TRUE)
t2_t0_GO_up <- read.csv("../../Data/results/GO/BP/t2_t0_GO_up.csv", row.names = 1, as.is = TRUE)

t3_t0_GO_down <- read.csv("../../Data/results/GO/BP/t3_t0_GO_down.csv", row.names = 1, as.is = TRUE)
t3_t0_GO_up <- read.csv("../../Data/results/GO/BP/t3_t0_GO_up.csv", row.names = 1, as.is = TRUE)



# DE genes
down_data_t1_t0 <- read.csv("../../Data/results/DE_genes/down_data_t1_t0.csv", row.names = 1)
up_data_t1_t0 <- read.csv("../../Data/results/DE_genes/up_data_t1_t0.csv", row.names = 1)

down_data_t2_t0 <- read.csv("../../Data/results/DE_genes/down_data_t2_t0.csv", row.names = 1)
up_data_t2_t0 <- read.csv("../../Data/results/DE_genes/up_data_t2_t0.csv", row.names = 1)

down_data_t3_t0 <- read.csv("../../Data/results/DE_genes/down_data_t3_t0.csv", row.names = 1)
up_data_t3_t0 <- read.csv("../../Data/results/DE_genes/up_data_t3_t0.csv", row.names = 1)
up_data_t3_t0 <- up_data_t3_t0[order(-up_data_t3_t0$logFC),]


amountOfTerms <- 6

# remove positive regulation of transcription fro... & negative regulation of transcription fro..
t3_t0_GO_up <- t3_t0_GO_up[-c(1,5),]

vector <- make_vector(t3_t0_GO_up[1:amountOfTerms,], up_data_t3_t0$geneName, HUGO)
genes_interest <- make_genes_interest(t3_t0_GO_up[1:amountOfTerms,], up_data_t3_t0$geneName, HUGO)

terms <- make_terms_dataframe(t3_t0_GO_up[1:amountOfTerms,],vector,"BP")

genes_interest
genes_test <- data.frame(ID = genes_interest)
genes_test$logFC <- up_data_t3_t0[match(genes_interest, up_data_t3_t0$geneName),"logFC"]




#-------------- make circle plot--------------------------

plot <- "top6_GO_terms.png"
#plot <- "top6_GO_terms_posNegRegulationREMOVED.png"
chord <- circle_plot(terms,genes_test,plot)
png(plot,height = 1000, width = 1500)
#GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 4,lfc.min = min(circ$logFC), lfc.max = max(circ$logFC))
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 4,lfc.min = min(chord[,ncol(chord)]), lfc.max = max(chord[,ncol(chord)]))
dev.off()







# load Swantje's 3 GO lists.
GO_list1 <- read.csv("GO_list1.csv", row.names = 1, as.is = TRUE)
GO_list2 <- read.csv("GO_list2.csv", row.names = 1, as.is = TRUE)
GO_list3 <- read.csv("GO_list3.csv", row.names = 1, as.is = TRUE)


### for GO_list1
vector <- make_vector(GO_list1, up_data_t3_t0$geneName, HUGO)
genes_interest <- make_genes_interest(GO_list1, up_data_t3_t0$geneName, HUGO)

terms <- make_terms_dataframe(GO_list1,vector,"BP")


genes_test <- data.frame(ID = genes_interest)
genes_test$logFC <- up_data_t3_t0[match(genes_interest, up_data_t3_t0$geneName),"logFC"]




#-------------- make circle plot--------------------------

plot <- "GO_list1.png"
#plot <- "top6_GO_terms_posNegRegulationREMOVED.png"
chord <- circle_plot(terms,genes_test,plot)
png(plot,height = 1000, width = 1500)
#GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 4,lfc.min = min(circ$logFC), lfc.max = max(circ$logFC))
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 4,lfc.min = min(chord[,ncol(chord)]), lfc.max = max(chord[,ncol(chord)]))
dev.off()





### for GO_list2
vector <- make_vector(GO_list2, up_data_t3_t0$geneName, HUGO)
genes_interest <- make_genes_interest(GO_list2, up_data_t3_t0$geneName, HUGO)

terms <- make_terms_dataframe(GO_list2,vector,"BP")


genes_test <- data.frame(ID = genes_interest)
genes_test$logFC <- up_data_t3_t0[match(genes_interest, up_data_t3_t0$geneName),"logFC"]




#-------------- make circle plot--------------------------

plot <- "GO_list2.png"
#plot <- "top6_GO_terms_posNegRegulationREMOVED.png"
chord <- circle_plot(terms,genes_test,plot)
png(plot,height = 1000, width = 1500)
#GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 4,lfc.min = min(circ$logFC), lfc.max = max(circ$logFC))
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 4,lfc.min = min(chord[,ncol(chord)]), lfc.max = max(chord[,ncol(chord)]))
dev.off()







### for GO_list3
vector <- make_vector(GO_list3, up_data_t3_t0$geneName, HUGO)
genes_interest <- make_genes_interest(GO_list3, up_data_t3_t0$geneName, HUGO)

terms <- make_terms_dataframe(GO_list3,vector,"BP")


genes_test <- data.frame(ID = genes_interest)
genes_test$logFC <- up_data_t3_t0[match(genes_interest, up_data_t3_t0$geneName),"logFC"]




#-------------- make circle plot--------------------------

plot <- "GO_list3.png"
#plot <- "top6_GO_terms_posNegRegulationREMOVED.png"
chord <- circle_plot(terms,genes_test,plot)
png(plot,height = 1000, width = 1500)
#GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 4,lfc.min = min(circ$logFC), lfc.max = max(circ$logFC))
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 4,lfc.min = min(chord[,ncol(chord)]), lfc.max = max(chord[,ncol(chord)]))
dev.off()

