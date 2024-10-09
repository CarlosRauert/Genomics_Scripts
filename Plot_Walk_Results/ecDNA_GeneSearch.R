library(gGnome)
library(genomation)

AdipocyteEnhancers <- readGeneric("Adipocyte.bed")

SearchEnhancersOnECDNA <- function(filepath){
  EcDNAs <- readRDS(filepath)

  Chromosomes=matrix("nn", ncol=15279, nrow=1, dimnames=list("Chromosome", seq(from=1,to=15279,by=1)))
  for (i in 1:15279){
    Chromosomes[1,i]<-as.character(seqnames(AdipocyteEnhancers[i]))
  }
  Observations=matrix(FALSE, ncol=15279, nrow=1, dimnames=list("Enhancer on ecDNA?", seq(from=1,to=15279,by=1)))
  ObservationsDF <- rbind(Observations, Chromosomes)
  
  
  for (i in 1:length(AdipocyteEnhancers)){
    Enhancer = AdipocyteEnhancers[i]
    Enhancer_ecDNAs <- EcDNAs$grl %&% Enhancer
    if (length(Enhancer_ecDNAs)>0){
      ObservationsDF[1,i]<-TRUE
    }
  }
  return(ObservationsDF)
}

SearchGenesOnECDNA <- function(filepath){
  EcDNAs <- readRDS(filepath)
  CGC <- readRDS("20190829cancergene_elements.rds")
  CGC_gr <- gr.sub(dt2gr(CGC[,c(1:3,6)]), 'chr', '')
  
  Chromosomes=matrix("nn", ncol=1059, nrow=1, dimnames=list("Chromosome", CGC_gr$geneName))
  for (i in 1:1059){
    Chromosomes[1,i]<-as.character(seqnames(CGC_gr[i]))
  }
  Observations=matrix(FALSE, ncol=1059, nrow=1, dimnames=list("Gene on ecDNA?", CGC_gr$geneName))
  ObservationsDF <- rbind(Observations,Chromosomes)
  
  
  for (i in 1:length(CGC_gr)){
    CGC_Gene = CGC_gr %Q% (geneName == CGC_gr$geneName[i])
    Gene_ecDNAs <- EcDNAs$grl %&% CGC_Gene
    if (length(Gene_ecDNAs)>0){
      ObservationsDF[1,i]<-TRUE
    }
  }
  return(ObservationsDF)
}  
A3LS_37<-SearchGenesOnECDNA("MDM2_Walks_Out/A3LSMDM2_Top3_walks_CNmin_37.25.rds")
A1L3_60 <- SearchGenesOnECDNA("MDM2_Walks_Out/A1L3MDM2_Top3_walks_CNmin_60.5.rds")
A1L3_Enhancers <- SearchEnhancersOnECDNA("MDM2_Walks_Out/A1L3MDM2_Top3_walks_CNmin_60.5.rds")
A3LS_Enhancers <- SearchEnhancersOnECDNA("MDM2_Walks_Out/A3LSMDM2_Top3_walks_CNmin_37.25.rds")

A3LS_37_all_Oncogenes <- SearchGenesOnECDNA(filepath = "MDM2_Walks_Out/A3LS/walks_circ_CNmin_37.25.rds")
A3LS_37_all_Enhancers <- SearchEnhancersOnECDNA(filepath = "MDM2_Walks_Out/A3LS/walks_circ_CNmin_37.25.rds")

