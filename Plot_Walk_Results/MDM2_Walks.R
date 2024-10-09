# MDM2 Walks

library(parallel)
library(gGnome)


# For all available Jabba Calls:

# Make List of Available Cases

CaseList <- c("A1KU_new","A1KW_new", "A1L0_new", "A1L2_new", "A1L3_new", "A2J4_new")

# Function that finds and plots circular amplicons

PlotecDNA <- function(CaseID){
  if(!dir.exists(paste0("MDM2_Walks_Out/",CaseID))){
    dir.create(paste0("MDM2_Walks_Out/",CaseID))
  }
  try({
    Ggraph <-  gG(jabba=paste0("Jabba_Calls_TCGA/",CaseID,"/jabba.simple.rds"))
    print(paste0("JabbA loaded for ", CaseID))
    gencode <- track.gencode(cached.dir="/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics", stack.gap = 1e5, cex.label = 0.8, height = 20)
    print("GENCODE loaded")
    CNmax <- max(na.omit(Ggraph$dt$cn))                                                                                                                                                     
    CNxt = CNmax *0.5                                                                                                                                                                       
    while (CNxt >=8) {                                                                                                                                                                      
      print(CNxt)
      highcopyX = Ggraph[cn>CNxt]        # subset to high copy only                                                                                                                                                   
      walks = highcopyX$walks()                                                                                                                                                           
      #saveRDS(walks, paste0('MDM2_Walks_Out/',CaseID,'/walks_CNmin_', CNxt,'.rds')) 
      print('Walks done')
      walks_circ <- walks[circular == T]            # ecDNA?                                                                                                                                      
      if (length(walks_circ) > 0) {
        print('Found a circular amplicon')
        # saveRDS(walks_circ, paste0('MDM2_Walks_Out/',CaseID,'/walksCirc_CNmin_', CNxt,'.rds'))                                                                                                 
        # we can choose the longest walk (most nodes traversed)                                                                                                                          
        walks_l = walks_circ[walk.id %in% names(sort(walks_circ$lengths, decreasing = T)[1:3])]                                                                                           
        pdf(file = paste0('MDM2_Walks_Out/',CaseID,'/walksCirc_CNmin_', CNxt,'.pdf'))                                                                                       
        plot(c(gencode, Ggraph$gt, walks_l$gtrack(name = "circular walks")), walks_l$footprint+2e5, cex.label=0.2)
        title(main=paste0("longest circular walks with CN greater than ",CNxt))
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

lapply(CaseList,PlotecDNA)