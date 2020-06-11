## ----pubmed1, fig.cap="Psychiatric disorders in PubMed. It has been obtained querying **psychiatric disorder [Title/Abstract] from 1955 to 2016**.", echo=FALSE, fig.wide = TRUE----

plot_pubmed_stats <-function(){

  library(ggplot2)
  
  data.file <- system.file(
    paste0("extdata", .Platform$file.sep, "psychiatricDisordersPubmed.csv"), 
    package="psygenet2r"
  )
  
  pmid <- read.delim(data.file, header=TRUE, sep=",")
  pmid <- pmid[pmid$year<2017 & pmid$year>1950,]
  
  pmid$year <- factor(pmid$year)
  labels <- as.integer(seq(1950, 2016, by=5))
  
  
  pdf("nPubPsychriatic.pdf",width=10,height=7)
  p <- ggplot(pmid, aes ( x = year, y = count ) ) +
    geom_bar ( stat = "identity", fill = "grey" ) +
    labs ( title = "Number of publications for psychiatric disorders in PubMed" , x = "year", y = "# of pmids") +
    theme_classic( ) + 
    scale_x_discrete(breaks=labels, labels=as.character(labels))+
    theme( plot.margin = grid::unit ( x = c ( 5, 15, 5, 15 ), units = "mm" ),
           axis.line = element_line ( size = 0.7, color = "black" ), 
           text = element_text ( size = 14 ) ,
           axis.text.x = element_text ( angle = 45, size = 11, hjust = 1 ) )
  
  p=grid.arrange(p,nrow = 1)
  dev.off()
  
}

use_psygenet2r <-function(genesOfInterest){
  
  library( psygenet2r )
  
  ## ----genes-----------------------------------------------------------------
  genesOfInterest <- c("ADCY2", "AKAP13", "ANK3", "ANKS1A", 
                       "ATP6V1G3", "ATXN1", "C11orf80", "C15orf53", "CACNA1C", 
                       "CACNA1D", "CACNB3", "CROT", "DLG2", "DNAJB4", "DUSP22", 
                       "FAM155A", "FLJ16124", "FSTL5", "GATA5", "GNA14", "GPR81", 
                       "HHAT", "IFI44", "ITIH3", "KDM5B", "KIF1A", "LOC150197", 
                       "MAD1L1", "MAPK10", "MCM9", "MSI2", "NFIX", "NGF", "NPAS3", 
                       "ODZ4", "PAPOLG", "PAX1", "PBRM1", "PTPRE", "PTPRT", 
                       "RASIP1", "RIMBP2", "RXRG", "SGCG", "SH3PXD2A", "SIPA1L2",
                       "SNX8", "SPERT", "STK39", "SYNE1", "THSD7A", "TNR", 
                       "TRANK1", "TRIM9", "UBE2E3", "UBR1", "ZMIZ1", "ZNF274")
  
  ## ----search_multiple-------------------------------------------------------
  m1 <- psygenetGene(
    gene     = genesOfInterest, 
    database = "ALL",
    verbose  = FALSE,
    warnings = FALSE
  )
  m1
  
  ## ----gene-disease, fig.height=8, fig.width=8, fig.cap = "Gene-Disease Association Network", fig.wide = TRUE----
  plot( m1 )
  
  ## ----gene-psy, fig.cap="Association type barplot according to psychiatric category", fig.wide = TRUE----
  geneAttrPlot( m1, type = "evidence index" )
  
  ## ----panther, fig.cap="Panther class analysis of the genes of interest.", message=FALSE, warning=FALSE, fig.wide = TRUE----
  pantherGraphic( genesOfInterest, "ALL")
  
  ## ----gene-disease-2, fig.cap="Gene-Disease Association Heatmap", fig.wide = TRUE----
  plot( m1, type="GDA heatmap")
  
  ## ----pubmed2, fig.cap="Publications that report each gene association with bipolar disorder", fig.wide = TRUE----
  plot( m1, name="bipolar disorder", type="publications")
  
  ## ----sentences1_query------------------------------------------------------
  m2 <- psygenetGeneSentences(
    geneList = genesOfInterest, 
    database = "ALL"
  )
  m2
  
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
  plot( xx )
  
  ## ----bpGenes, fig.cap="Barplot: Genes associated to each of the psychiatric disorders", fig.wide = TRUE----
  geneAttrPlot( m1, type = "disease category" )
  
  ## ----bpDis, fig.cap="Barplot: CUIs and psychiatric categories associated to each gene", fig.wide = TRUE----
  geneAttrPlot( m1, type = "gene" )

  
}