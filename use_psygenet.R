## ----pubmed1, fig.cap="Psychiatric disorders in PubMed. 
## It has been obtained querying **psychiatric disorder [Title/Abstract] from 1955 to 2016**.", 

use_psygenet2r_disease <-function(genesOfInterest,id){
  
  library( psygenet2r )
  
  ## ----genes-----------------------------------------------------------------
  #genesOfInterest <- c("ADCY2", "AKAP13", "ANK3", "ANKS1A", 
  #                     "ATP6V1G3", "ATXN1", "C11orf80", "C15orf53", "CACNA1C", 
  #                     "CACNA1D", "CACNB3", "CROT", "DLG2", "DNAJB4", "DUSP22", 
  #                     "FAM155A", "FLJ16124", "FSTL5", "GATA5", "GNA14", "GPR81", 
  #                     "HHAT", "IFI44", "ITIH3", "KDM5B", "KIF1A", "LOC150197", 
  #                     "MAD1L1", "MAPK10", "MCM9", "MSI2", "NFIX", "NGF", "NPAS3", 
  #                     "ODZ4", "PAPOLG", "PAX1", "PBRM1", "PTPRE", "PTPRT", 
  #                     "RASIP1", "RIMBP2", "RXRG", "SGCG", "SH3PXD2A", "SIPA1L2",
  #                     "SNX8", "SPERT", "STK39", "SYNE1", "THSD7A", "TNR", 
  #                     "TRANK1", "TRIM9", "UBE2E3", "UBR1", "ZMIZ1", "ZNF274")
  
  ## ----search_multiple-------------------------------------------------------
  m1 <- psygenetGene(
    gene     = genesOfInterest, 
    database = "ALL",
    verbose  = FALSE,
    warnings = FALSE
  )
  #m1
  
  ## ----gene-disease, fig.height=8, fig.width=8, fig.cap = "Gene-Disease Association Network", fig.wide = TRUE----
  fName = sprintf("./gene_disease_%d.pdf",id)
  pdf(fName,width=10,height=8)
  plot( m1 )
  dev.off()
  
  ## ----gene-psy, fig.cap="Association type barplot according to psychiatric category", fig.wide = TRUE----
  pdf("./gene_psy.pdf",width=10,height=8)
  geneAttrPlot( m1, type = "evidence index" )
  dev.off()
  
  ## ----panther, fig.cap="Panther class analysis of the genes of interest.", message=FALSE, warning=FALSE, fig.wide = TRUE----
  #pdf("./panther.pdf",width=10,height=8)
  #pantherGraphic( genesOfInterest, "ALL")
  #dev.off()
  
  ## ----gene-disease-2, fig.cap="Gene-Disease Association Heatmap", fig.wide = TRUE----
  pdf("./gene_disease_2.pdf",width=10,height=8)
  plot( m1, type="GDA heatmap")
  dev.off()
  
  ## ----pubmed2, fig.cap="Publications that report each gene association with bipolar disorder", fig.wide = TRUE----
  pdf("./gene_bipolar.pdf",width=10,height=8)
  plot( m1, name="bipolar disorder", type="publications")
  dev.off()
  
}  
  
use_psygenet2r_sentence <-function(genesOfInterest){
  
  ## ----sentences1_query------------------------------------------------------
  m2 <- psygenetGeneSentences(
    geneList = genesOfInterest, 
    database = "ALL"
  )
  #m2
  
  ## ----sentences2_extraction, warnings=FALSE---------------------------------
  sentences <- extractSentences( m2, 
                                 disorder = "bipolar disorder" )
  head(sentences$PUBMED_ID)
  
  ## ----jaccard_1, warning=FALSE----------------------------------------------
  xx <- jaccardEstimation( genesOfInterest, "bipolar disorder", 
                           database = "ALL", nboot = 500 )
  xx
  
  ## ----jaccard_2-------------------------------------------------------------
  extract( xx )
  
  ## ----jaccard_3, warning=FALSE----------------------------------------------
  xx <- jaccardEstimation( genesOfInterest, 
                           database = "ALL", nboot = 500 )
  
  ## ----jacc, fig.cap="Bar-plot where the Jaccard Index of each comparison between the list of genes of interest and PsyGeNET's diseases is shown.", fig.wide = TRUE----
  pdf("./jacc.pdf",width=10,height=8)
  plot( xx )
  dev.off()
  
  ## ----bpGenes, fig.cap="Barplot: Genes associated to each of the psychiatric disorders", fig.wide = TRUE----
  pdf("./bpGenes.pdf",width=10,height=8)
  geneAttrPlot( m1, type = "disease category" )
  dev.off()
  
  ## ----bpDis, fig.cap="Barplot: CUIs and psychiatric categories associated to each gene", fig.wide = TRUE----
  pdf("./bpDis.pdf",width=10,height=8)
  geneAttrPlot( m1, type = "gene" )
  dev.off()
  
}