library(gGnome)
library(parallel)
gt.ge <- readRDS('/inputdata/20240626jabba_walks/gt_ge.rds')
CGC <- readRDS('/inputdata/20190829cancergene_elements.rds')
CGC_gr <- gr.sub(dt2gr(CGC[,c(1:3,6)]), 'chr', '')
#xT = 'Lu'
mclapply(c('BRCA', 'NSCLC'), function(xT){
  # xG = 'CCND1'
  mclapply(c('CCND1', 'ERBB2', 'MYC', 'FGFR1'), function(xG){
    print(xG) 
    CGC_xGgr = CGC_gr %Q% (geneName == xG)  
    CNlist_dt <- readRDS(paste0('/inputdata/20240626jabba_walks/in/20240107_', xT, '_', xG, '_CNlist_dt.rds'))   # list of files to run on                                                                                                                           
    fls <- unique(unique(CNlist_dt$Tumor_Sample_Barcode))
    #fx = fls[1]
    CNlist <- mclapply(fls, function(fx){                                                                                                                                                     
      print(fx)                                                                                                                                                                               
      jabbaX <- gG(jabba = paste0('/inputdata/chooetal2023jabbacalls/20240411jabba/', fx, '/jabba.simple.rds'))        # read in jabba                                                                                                                                        
      CNmax <- max(na.omit(jabbaX$dt$cn))                                                                                                                                                     
      CNxt = CNmax *0.5                                                                                                                                                                       
      while (CNxt >=8) {                                                                                                                                                                      
        print(CNxt)
        highcopyX = jabbaX[cn>CNxt]        # subst to high copy only                                                                                                                                                   
        try(expr ={                                                                                                                                                                           
          walks = highcopyX$walks()                                                                                                                                                           
          saveRDS(walks, paste0('/inputdata/20240626jabba_walks/out/', xT,'_', xG, '_', fx, 'walks_cnLrg', CNxt,'.rds')) 
          print('Walks done')
          walks_circ <- walks[circular == T]            # ecDNA?                                                                                                                                      
          if (length(walks_circ) > 0) {
            print('Found a circular amplicon')
            saveRDS(walks_circ, paste0('/inputdata/20240626jabba_walks/out/', xT,'_', xG, '_', fx, 'walksCirc_cnLrg', CNxt,'.rds'))                                                                                                 
            walks_circs_dt <- gr2dt(walks_circ$grl)                                                                                                                                           
            walks_circs_dtgr <- dt2gr(walks_circs_dt)                                                                                                                                         
            # get all with xG                                                                                                                                                                 
            walks_circs_dtgr_xGrdt <- gr2dt(walks_circs_dtgr %&% CGC_xGgr) # MDM2 overlap                                                                                                             
            walks_circs_dtgr_xGrdtan <- gr2dt(walks_circs_dtgr %$% CGC_gr) #                                                                                                                  
            unique(walks_circs_dtgr_xGrdtan$geneName)                                                                                                                                         
            # walks_circ_sdt <- walks_circ_s$dt                                                                                                                                               
            ## we can choose the longest walk (most nodes traversed)                                                                                                                          
            walks_l = walks_circ[walk.id %in% names(sort(walks_circ$lengths, decreasing = T)[1:3])]                                                                                           
            walks_circ_xGr <- walks_circ %&% CGC_xGgr                                                                                                                                         
            if (length(walks_circ_xGr)>0) {                                                                                                                                                   
              saveRDS(walks_circ_xGr, paste0('/inputdata/20240626jabba_walks/out/', xT,'_',xG, '_', fx, 'OncogeneWalksCirc_cnLrg', CNxt,'.rds'))                                                                                   
              pdf(file = paste0('/inputdata/20240626jabba_walks/out/', xT,'_',xG, '_', fx, 'walksCirc_cnLrg', CNxt,'.pdf'), paper = 'legal')                                                                                       
              plot(c(gt.ge, jabbaX$gt, walks_l$gtrack(name = "Longest walks")), walks_l$footprint+2e5, cex.label = 0.15)                                                                      
              plot(c(gt.ge, jabbaX$gt, walks_circ_xGr$gtrack(name = "Oncogene walks", height = 40)), walks_circ_xGr$footprint+2e5, cex.label = 0.1)                                           
              dev.off()                                                                                                                                                                       
              break                                                                                                                                                                           
            } else {                                                                                                                                                                          
              pdf(file = paste0('/inputdata/20240626jabba_walks/out/', xT,'_', xG, '_',fx, 'walksCirc_cnLrg', CNxt,'.pdf'), paper = 'legal')                                                                                       
              plot(c(gt.ge, jabbaX$gt, walks_l$gtrack(name = "Longest walks")), walks_l$footprint+2e5, cex.label = 0.15)                                                                      
              dev.off()                                                                                                                                                                       
            }                                                                                                                                                                                 
          }                                                                                                                                                                                   
          CNxt = CNxt *2/3                                                                                                                                                                    
          if (CNxt <7) {                                                                                                                                                                      
            break                                                                                                                                                                             
          }                                                                                                                                                                                   
          if (CNxt < 10 & CNxt > 7.5) {                                                                                                                                                       
            CNxt=8                                                                                                                                                                            
          }                                                                                                                                                                                   
        })                                                                                                                                                                                    
      }
    }, mc.cores = 3, mc.preschedule = F)
  }, mc.cores = 4, mc.preschedule = F,mc.allow.recursive = TRUE) 
}, mc.cores = 2, mc.preschedule = F,mc.allow.recursive = TRUE)

