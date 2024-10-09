library(gGnome)

AllEcDNA <- function(walks){
  All <- readRDS(walks)
  walks_circ <- All[circular == T] # ecDNA?
}

All_A1L3 <- AllEcDNA("MDM2_Walks_Out/A1L3/walks_CNmin_60.5.rds")

saveRDS(All_A1L3,"MDM2_Walks_Out/A1L3/walks_circ_CNmin_60.5.rds")

gencode <- track.gencode(cached.dir="/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics", stack.gap = 1e5, cex.label = 20, height = 20)
gencode$legend=FALSE
plot(c(gencode, All_A1L3$gtrack(name = "circular walks")), All_A1L3$footprint+2e5, cex.label=0.8)

All_A3LS_CN37 <- AllEcDNA("MDM2_Walks_Out/A3LS/walks_CNmin_37.25.rds")

saveRDS(All_A3LS_CN37,"MDM2_Walks_Out/A3LS/walks_circ_CNmin_37.25.rds")


plot(c(gencode, All_A3LS_CN37$gtrack(name = "circular walks")), All_A3LS_CN37$footprint+2e5, cex.label=0.8)

