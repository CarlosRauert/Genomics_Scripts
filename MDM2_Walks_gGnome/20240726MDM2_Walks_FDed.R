# MDM2 Walks

library(parallel)
library(gGnome)

CGC <- readRDS('20190829cancergene_elements.rds') # cancer gene census genes https://cancer.sanger.ac.uk/census
CGC_gr <- gr.sub(dt2gr(CGC[,c(1:3,6)]), 'chr', '')
xG = 'MDM2'
CGC_xGgr = CGC_gr %Q% (geneName == xG)  
MDM2_gr <- CGC_xGgr
# now use   walks_circ_xGr <- walks_circ %&% CGC_xGgr  (also see https://github.com/mskilab-org/gUtils)


# For all available Jabba Calls:

# Make List of Available Cases

# CaseList <- read_tsv("Jabba_Calls_TCGA/JabBa_TCGA_Liposarc.tsv")$pair

# Function that finds and plots circular amplicons

PlotecDNA <- function(CaseID){
  if(!dir.exists(paste0("/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics/MDM2_Walks_Out/",CaseID))){
    dir.create(paste0("/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics/MDM2_Walks_Out/",CaseID))
  }
  try({
    Ggraph <-  gG(jabba=paste0("/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics/Jabba_Calls_TCGA/",CaseID,"/jabba.simple.rds"))
    print(paste0("JabbA loaded for ", CaseID))
    CNmax <- max(na.omit(Ggraph$dt$cn))                                                                                                                                                     
    CNxt = CNmax *0.5                                                                                                                                                                       
    while (CNxt >=8) {                                                                                                                                                                      
      print(CNxt)
      highcopyX = Ggraph[cn>CNxt]        # subset to high copy only                                                                                                                                                   
      walks = highcopyX[1:length(highcopyX)]$walks(verbose=TRUE)                                                                                                                                                           
      saveRDS(walks, paste0('/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics/MDM2_Walks_Out/',CaseID,'/walks_CNmin_', CNxt,'.rds')) 
      print('Walks done')
      walks_circ <- walks[circular == T] # ecDNA?
      print("checked for circular walks")
      walks_circ_MDM2_gr <- walks_circ$grl %&% MDM2_gr # contain MDM2?
      walks_circ_MDM2 <- walks_circ[as.integer(names(walks_circ_MDM2_gr))]
      if (length(walks_circ_MDM2) > 0) {
        print('Found circular MDM2 amplicon')
        # we can choose the longest walk (most nodes traversed)                                                                                                                          
        walks_l = walks_circ_MDM2[walk.id %in% names(sort(walks_circ_MDM2$lengths, decreasing = T)[1:3])]
        saveRDS(walks_l, paste0("/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics/MDM2_Walks_Out/",CaseID,"MDM2_Top3_walks_CNmin_", CNxt,".rds"))
        pdf(file = paste0('/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics/MDM2_Walks_Out/',CaseID,'/walksCircMDM2_CNmin_', CNxt,'.pdf'))                                                                                       
        plot(c(gencode, Ggraph$gt, walks_l$gtrack(name = "circular walks")), walks_l$footprint+2e5, cex.label=0.2)
        title(main=paste0("MDM2 circular walks with CN greater than ",CNxt))
        par(cex.lab=0.1)
        dev.off()                                                                                                                                                                       
      }
      CNxt <- CNxt/2                                                                                                                                                                   
      if (CNxt <7) {                                                                                                                                                                      
        break                                                                                                                                                                             
      }
      if (CNxt < 10 & CNxt > 7.5) {                                                                                                                                                       
         CNxt=8                                                                                                                                                                            
       }  
    }     
  })
}
# Run
gencode <- track.gencode(cached.dir="/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics", stack.gap = 1e5, cex.label = 0.8, height = 20)

PlotecDNA("A1KU")

lapply(CaseList,PlotecDNA)