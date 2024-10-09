library(readr)
library(genomation)
library(gGnome)
library(pbmcapply)

CaseTable <- read_tsv("Jabba_Calls_TCGA/JabBa_TCGA_Liposarc.tsv")
CaseList <- CaseTable$pair

SearchGenesOnECDNA <- function(filepath){
  EcDNAs <- readRDS(filepath)
  walks_circ <- EcDNAs[circular == T] # ecDNA?
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

ReadAllEcDNAs <- function(CaseID){
  FileList <- list.files(paste0("MDM2_Walks_Out/",CaseID))
  FileList <- grep("walks_CNmin_",FileList, value=TRUE)
  FileList <- paste0("MDM2_Walks_Out/",CaseID,"/",FileList)
  FileList
  OncogenesInEcDNA_Sample <- list()
  for (i in 1:length(FileList)){
    print(paste0("Looking through ecDNAs in File ", FileList[i]))
    Observations_i <- SearchGenesOnECDNA(FileList[i])
    GenesOnEcDNA <- names(which(Observations_i[1,]==T))
    OncogenesInEcDNA_Sample <- append(OncogenesInEcDNA_Sample, GenesOnEcDNA)
    OncogenesInEcDNA_Sample <- unique(OncogenesInEcDNA_Sample)
  }
  return(OncogenesInEcDNA_Sample)
}

MakeEcDNA_Table <- function(CaseList){
  GenesList<-pbmclapply(CaseList, ReadAllEcDNAs, mc.cores=14)
  EcDNA_Table <- rbind(CaseList, GenesList)
  return(EcDNA_Table)
}  

Table <- MakeEcDNA_Table(CaseList)
