# A1KU test with modified functions

# Walks Maker thats ready to be passed to shell script, starts from specified CN

Args <- commandArgs(trailingOnly=TRUE)
CaseID <- Args[1]
CNxt <- as.integer(Args[2])

library(parallel)
library(gGnome)
print(CaseID)
print(CNxt)
print("test")
CGC <- readRDS('/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics/20190829cancergene_elements.rds') # cancer gene census genes https://cancer.sanger.ac.uk/census
CGC_gr <- gr.sub(dt2gr(CGC[,c(1:3,6)]), 'chr', '')
xG = 'MDM2'
CGC_xGgr = CGC_gr %Q% (geneName == xG)  
MDM2_gr <- CGC_xGgr

gencode <- track.gencode(cached.dir="/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics", stack.gap = 1e5, cex.label = 0.8, height = 20)

sparse_subset = function(A, B, strict = FALSE, chunksize = 100, quiet = FALSE)
{
  nz = Matrix::colSums(as.matrix(A)!=0, 1)>0
  
  if (is.null(dim(A)) | is.null(dim(B)))
    return(NULL)
  
  C = Matrix::sparseMatrix(i = c(), j = c(), dims = c(nrow(A), nrow(B)))
  
  for (i in seq(1, nrow(A), chunksize))
  {
    ixA = i:min(nrow(A), i+chunksize-1)
    for (j in seq(1, nrow(B), chunksize))
    {
      ixB = j:min(nrow(B), j+chunksize-1)
      
      if (length(ixA)>0 & length(ixB)>0 & !quiet)
        cat(sprintf('\t interval A %s to %s (%d) \t interval B %d to %d (%d)\n', ixA[1], ixA[length(ixA)], nrow(A), ixB[1], ixB[length(ixB)], nrow(B)))
      if (strict)
        C[ixA, ixB] = (sign((A[ixA, , drop = FALSE]!=0)) %*% sign(t(B[ixB, , drop = FALSE]!=0))) * (sign((A[ixA, , drop = FALSE]==0)) %*% sign(t(B[ixB, , drop = FALSE]!=0))>0)
      else
        C[ixA, ixB] = (sign(A[ixA, nz, drop = FALSE]!=0) %*% sign(t(B[ixB, nz, drop = FALSE]==0)))==0
    }
  }
  
  return(C)
}

convex.basis2 = function(A,
                         interval = 80,
                         chunksize = 100,
                         exclude.basis = NULL,
                         exclude.range = NULL,
                         maxchunks = Inf,
                         verbose = F){
  ZERO = 1e-8;
  remaining = 1:nrow(A);
  iter = 0;
  i = 0;
  #    order = c()
  numelmos = c()
  K_i = I = as(diag(rep(1, ncol(A))), 'sparseMatrix');
  #    A_i = as(A %*% K_i, 'sparseMatrix');
  K_i = I = diag(rep(1, ncol(A)))
  A_i = A %*% K_i
  print("Convex_Basis is edited")
  if (!is.null(exclude.basis)){
    exclude.basis = sign(exclude.basis)
    exclude.basis = exclude.basis[rowSums(exclude.basis)>0, ]
    if (nrow(exclude.basis) == 0){
      exclude.basis = NULL
    }
  }
  
  if (!is.null(exclude.range)){
    exclude.range = sign(exclude.range)
    exclude.range = exclude.range[rowSums(exclude.range)>0, ]
    if (nrow(exclude.range) == 0){
      exclude.range = NULL
    }
  }
  
  # vector to help rescale matrix (avoid numerical issues)
  mp  = apply(abs(A), 1, min); # minimum value of each column
  mp[mp[ZERO]] = 1; # columns with zero minimum get scale "1"
  
  st = Sys.time()
  # iterate through rows of A, "canceling" them out
  while (length(remaining)>0){   
    ## TODO figure out why we have to check this so many times
    if (nrow(K_i)==0 | ncol(K_i)==0){
      return(matrix())
    }
    
    iter = iter+1;
    K_last = K_i;
    
    if (verbose){
      print(Sys.time() - st)
    }
    
    if (verbose){
      cat('Iter ', iter, '(of',  nrow(A_i),  ') Num basis vectors: ', nrow(K_i), " Num active components: ", sum(Matrix::rowSums(K_i!=0)), "\n")
    }
    
    i = remaining[which.min(Matrix::rowSums(A_i[remaining,, drop = FALSE]>=ZERO)*Matrix::rowSums(A_i[remaining,, drop = FALSE]<=(-ZERO)))]  # chose "cheapest" rows
    
    remaining = setdiff(remaining, i);
    #        order = c(order, i);
    
    zero_elements = which(abs(A_i[i, ]) <= ZERO);
    K_i1 = K_last[zero_elements, , drop = FALSE];  ## K_i1 = rows of K_last that are already orthogonal to row i of A
    K_i2 = NULL; ## K_i1 = will store positive combs of K_last rows that are orthogonal to row i of A (will compute these below)
    
    pos_elements = which(A_i[i, ]>ZERO)
    neg_elements = which(A_i[i, ]<(-ZERO))
    
    if (verbose){
      cat('Iter ', iter, " Row ", i, ":", length(zero_elements), " zero elements ", length(pos_elements), " pos elements ", length(neg_elements), " neg elements \n")
    }
    
    if (length(pos_elements)>0 & length(neg_elements)>0)
      for (m in seq(1, length(pos_elements), interval))
        for (l in seq(1, length(neg_elements), interval)){
          ind_pos = c(m:min(c(m+interval, length(pos_elements))))
          ind_neg = c(l:min(c(l+interval, length(neg_elements))))
          
          indpairs = cbind(rep(pos_elements[ind_pos], length(ind_neg)),
                           rep(neg_elements[ind_neg], each = length(ind_pos))); # cartesian product of ind_pos and ind_neg
          pix = rep(1:nrow(indpairs), 2)
          ix = c(indpairs[,1], indpairs[,2])
          #                coeff = c(-A_i[i, indpairs[,2]], A_i[i, indpairs[,1]])  ## dealing with Matrix ghost
          coeff = c(-A_i[i, ][indpairs[,2]], A_i[i, ][indpairs[,1]])  ##
          combs = Matrix::sparseMatrix(pix, ix, x = coeff, dims = c(nrow(indpairs), nrow(K_last)))
          combs[cbind(pix, ix)] = coeff;
          
          H = combs %*% K_last;
          
          #remove duplicated rows in H (with respect to sparsity)
          #H = H[!duplicated(as.matrix(H)>ZERO), , drop = F];
          
          # remove rows in H that have subsets in H (with respect to sparsity) ..
          if ((as.numeric(nrow(H))*as.numeric(nrow(H)))>maxchunks){
            print('Exceeding maximum number of chunks in convex.basis computation')
            stop('Exceeding maximum number of chunks in convex.basis computation')
          }
          if(nrow(abs(H)>ZERO)>0){
            keep = which(Matrix::colSums(sparse_subset(abs(H)>ZERO, abs(H)>ZERO, chunksize = chunksize, quiet = !verbose))<=1) # <=1 since every H is its own subset
            H = H[keep, , drop = FALSE]
          }
          # remove rows in H that have subsets in K_i2
          if (!is.null(K_i2))
            if (nrow(K_i2)>0){
              if ((as.numeric(nrow(K_i2))*as.numeric(nrow(H)))>maxchunks){
                print('Exceeding maximum number of chunks in convex.basis computation')
                stop('Exceeding maximum number of chunks in convex.basis computation')
              }
              if(nrow(abs(H)>ZERO)>0){
                keep = which(Matrix::colSums(sparse_subset(abs(K_i2)>ZERO, abs(H)>ZERO, chunksize = chunksize, quiet = !verbose))==0)
                H = H[keep, , drop = FALSE]
              }
            }
          
          #remove rows in H that have subsets in K_i1
          if (!is.null(K_i1))
            if (nrow(K_i1)>0){
              if ((as.numeric(nrow(K_i1))*as.numeric(nrow(H)))>maxchunks)
              {
                print('Exceeding maximum number of chunks in convex.basis computation')
                stop('Exceeding maximum number of chunks in convex.basis computation')
              }
              if(nrow(abs(H)>ZERO)>0){
                keep = which(Matrix::colSums(sparse_subset(abs(K_i1)>ZERO, abs(H)>ZERO, chunksize = chunksize, quiet = !verbose))==0)
                H = H[keep, , drop = FALSE]
              }
            }
          
          # maintain numerical stability
          if ((iter %% 10)==0){
            H = diag(1/apply(abs(H), 1, max)) %*% H
            
            #                K_i2 = rBind(K_i2, H)
          }
          K_i2 = rbind(K_i2, as.matrix(H))
        }
    
    #        K_i = rBind(K_i1, K_i2)
    K_i = rbind(K_i1, K_i2) ## new basis set
    
    if (nrow(K_i)==0){
      return(matrix())
    }
    
    ## only keep vectors that fail to intersect all vectors "exclude" in matrix
    if (!is.null(exclude.basis)) {
      if ((as.numeric(nrow(exclude.basis))*as.numeric(nrow(K_i)))>maxchunks){
        print('Exceeding maximum number of chunks in convex.basis computation')
        stop('Exceeding maximum number of chunks in convex.basis computation')
      }
      keep = Matrix::colSums(sparse_subset(exclude.basis>0, K_i>ZERO))==0
      if (verbose){
        cat('Applying basis exclusion and removing', sum(keep==0), 'basis vectors\n')
      }
      K_i = K_i[keep, , drop = F]
    }
    
    ## only keep vectors that fail to intersect all vectors "exclude" in matrix
    if (!is.null(exclude.range)){
      A_i_abs = abs(A) %*% t(K_i)
      if ((as.numeric(nrow(exclude.range))*as.numeric*ncol(A_i_abs))>maxchunks){
        print('Exceeding maximum number of chunks in convex.basis computation')
        stop('Exceeding maximum number of chunks in convex.basis computation')
      }
      keep = Matrix::colSums(sparse_subset(exclude.range>0, t(A_i_abs), quiet = !verbose))==0
      if (verbose){
        cat('Applying range exclusion and removing', sum(keep==0), 'basis vectors\n')
      }
      K_i = K_i[keep, , drop = F]
    }
    
    A_i = A %*% t(K_i)
  }
  
  return(t(K_i))
}

all.paths2 = function(A, all = F, ALL = F, sources = c(), sinks = c(), source.vertices = sources, sink.vertices = sinks,
                     exclude = NULL, ## specifies illegal subpaths, all such paths / cycles and
                     ## their supersets will be excluded, specified as k x nrow(A) matrix of vertex sets
                     verbose = FALSE,...)
{
  blank.vertices = Matrix::which(Matrix::rowSums(A)==0 & Matrix::colSums(A)==0)
  
  if (ALL)
    all = T
  
  if (all)
  {
    source.vertices = Matrix::which(Matrix::rowSums(A)>0 & Matrix::colSums(A)==0)
    sink.vertices = Matrix::which(Matrix::colSums(A)>0 & Matrix::rowSums(A)==0)
  }
  
  out = list(cycles = NULL, paths = NULL)
  
  node.ix = which(Matrix::rowSums(A!=0)>0 | Matrix::colSums(A!=0)>0)
  if (length(node.ix)==0)
    return(out)
  
  A = A[node.ix, node.ix]
  
  if (!is.null(exclude))
    exclude = sign(abs(exclude[, node.ix]))
  
  ij = Matrix::which(A!=0, arr.ind = T)
  B = Matrix::sparseMatrix(c(ij[,1], ij[,2]), rep(1:nrow(ij), 2), x = rep(c(-1, 1), each = nrow(ij)), dims = c(nrow(A), nrow(ij)))
  I = diag(rep(1, nrow(A)))
  
  source.vertices = setdiff(match(source.vertices, node.ix), NA)
  sink.vertices = setdiff(match(sink.vertices, node.ix), NA)
  
  B2 = Matrix::cbind2(Matrix::cbind2(B, I[, source.vertices, drop = FALSE]), -I[, sink.vertices, drop = FALSE])
  
  if (verbose)
    cat(sprintf('Computing paths for %s vertices and %s edges\n', nrow(B2), ncol(B2)))
  
  K = convex.basis2(B2, verbose = verbose, exclude.range = exclude, ...)
  
  if (all(is.na(K)))
    return(out)
  
  K = K[, Matrix::colSums(K[1:ncol(B), ,drop = FALSE])!=0, drop = FALSE] ## remove any pure source to sink paths
  
  is.cyc = Matrix::colSums(B %*% K[1:ncol(B), ,drop = FALSE]!=0)==0
  
  
  out$cycles = lapply(which(is.cyc),
                      function(i)
                      {
                        k = which(K[1:ncol(B), i]!=0)
                        v.all = unique(as.vector(ij[k, , drop = FALSE]))
                        sG = graph.edgelist(ij[k, , drop = FALSE])
                        tmp.v = v.all[c(1,length(v.all))]
                        p.fwd = get.shortest.paths(sG, tmp.v[1], tmp.v[2])
                        p.bwd = get.shortest.paths(sG, tmp.v[2], tmp.v[1])
                        return(node.ix[unique(unlist(c(p.fwd, p.bwd)))])
                      })
  
  out$paths = lapply(which(!is.cyc),
                     function(i)
                     {
                       k = K[1:ncol(B), i]
                       eix = which(k!=0)
                       v.all = unique(as.vector(ij[eix, , drop = FALSE]))
                       sG = graph.edgelist(ij[eix, , drop = FALSE])
                       io = B %*% k
                       v.in = Matrix::which(io<0)[1]
                       v.out = Matrix::which(io>0)[1]
                       return(node.ix[unlist(get.shortest.paths(sG, v.in, v.out))])
                     })
  
  if (length(out$cycles)>0)
  {
    tmp.cix = cbind(unlist(lapply(1:length(out$cycles), function(x) rep(x, length(out$cycles[[x]])))), unlist(out$cycles))
    out$cycles = out$cycles[!duplicated(as.matrix(Matrix::sparseMatrix(tmp.cix[,1], tmp.cix[,2], x = 1)))]
  }
  
  if (length(out$paths)>0)
  {
    tmp.pix = cbind(unlist(lapply(1:length(out$paths), function(x) rep(x, length(out$paths[[x]])))), unlist(out$paths))
    out$paths = out$paths[!duplicated(as.matrix(Matrix::sparseMatrix(tmp.pix[,1], tmp.pix[,2], x = 1)))]
  }
  
  if (ALL & length(blank.vertices)>0)
    out$paths = c(out$paths, lapply(blank.vertices, identity))
  
  return(out)
}
gGraph$set("public", "walks2", function(field = NULL, greedy = FALSE, verbose = FALSE) {
  
  if (greedy==TRUE){
    stop('Greedy walks not yet implemented - please stay tuned')
  }
  A = self$adj
  colnames(A) = rownames(A) = self$gr$snode.id
  
  
  lleft = which(self$nodes$loose.left)
  lright = which(self$nodes$loose.right)
  
  llefti = self$queryLookup(lleft)$index
  lrighti = self$queryLookup(lright)$index
  
  lleftir = self$queryLookup(-lleft)$index
  lrightir = self$queryLookup(-lright)$index
  
  sources = c(llefti, lrightir)
  sinks = c(lrighti, lleftir)
  
  ap = all.paths2(A, sources = sources, sinks = sinks, verbose = verbose)
  
  ## first make paths from any nodes that are both
  ## sources and sinks
  tmp = self$gr$snode.id[intersect(sources, sinks)]
  paths = split(tmp, seq_along(tmp))
  
  paths = c(paths, lapply(ap$paths, function(x) self$gr$snode.id[x]))
  cycles = lapply(ap$cycles, function(x) self$gr$snode.id[x])
  
  ## dedup reciprocals
  ## gets tricky for cycles since may be
  ## out of phase, however the (unordered) set of vertices / antivertices
  ## will be unique to that path / cycle and it's reciprocal
  ## so we don't need to worry about the sequence to match these up
  ## can just sort and match
  
  
  ## label each path / cycle
  pstr = sapply(paths, function(x) paste(sort(x), collapse = ' '))
  cstr = sapply(cycles, function(x) paste(sort(x), collapse = ' '))
  
  ## label reverse complement version
  rpstr = sapply(paths, function(x) paste(sort(-x), collapse = ' '))
  rcstr = sapply(cycles, function(x) paste(sort(-x), collapse = ' '))
  
  
  ## match up paths and their reverse complement
  ped = cbind(seq_along(rpstr), match(pstr, rpstr))
  pcl = igraph::clusters(igraph::graph.edgelist(ped), 'weak')$membership
  
  ## match up cycles and their reverse complement
  ## use vgrep here since cycles generated by all paths
  ## may be "out of phase" with their reverse complements
  ced = cbind(seq_along(rcstr), match(cstr, rcstr))
  ccl = igraph::clusters(igraph::graph.edgelist(ced), 'weak')$membership
  
  paths = paths[!duplicated(pcl)]
  cycles = cycles[!duplicated(ccl)]
  
  circular = c(rep(FALSE, length(paths)),
               rep(TRUE, length(cycles)))
  
  return(gWalk$new(snode.id = c(paths, cycles), graph = self, circular = circular))
})

PlotecDNA <- function(CaseID){
  if(!dir.exists(paste0("/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics/MDM2_Walks_Out/",CaseID))){
    dir.create(paste0("/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics/MDM2_Walks_Out/",CaseID))
  }
  try({
    Ggraph <-  gG(jabba=paste0("/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics/Jabba_Calls_TCGA/",CaseID,"/jabba.simple.rds"))
    print(paste0("JabbA loaded for ", CaseID))
    if(CNxt==999){
      CNmax <- max(na.omit(Ggraph$dt$cn))                                                                                                                                                     
      CNxt = CNmax *0.5
      print("Starting from the top")
    }
    print(paste0("Starting with ", CNxt))                                                                                                                                                                 
    while (CNxt >=8) {                                                                                                                                                                      
      print(CNxt)
      highcopyX = Ggraph[cn>CNxt]        # subset to high copy only                                                                                                                                                   
      walks = highcopyX[1:length(highcopyX)]$walks2(verbose=TRUE)                                                                                                                                                           
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
        plot(c(gencode, Ggraph$gt, walks_l$gtrack(name = "circular walks")), walks_l$footprint+2e5, cex.label=0.8)
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
PlotecDNA(CaseID)

print("done")
