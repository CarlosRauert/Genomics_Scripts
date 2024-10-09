library(gGnome)

# Functions for testing:

all.paths = function(A, all = F, ALL = F, sources = c(), sinks = c(), source.vertices = sources, sink.vertices = sinks,
                     exclude = NULL, ## specifies illegal subpaths, all such paths / cycles and
                     ## their supersets will be excluded, specified as k x nrow(A) matrix of vertex sets
                     verbose = TRUE,...)
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
  
  K = convex.basis(B2, verbose = verbose, exclude.range = exclude, ...)
  
  if (all(is.na(K)))
    return(out)
  
  K = K[, Matrix::colSums(K[1:ncol(B), ,drop = FALSE])!=0, drop = FALSE] ## remove any pure source to sink paths
  
  is.cyc = Matrix::colSums(B %*% K[1:ncol(B), ,drop = FALSE]!=0)==0
  
  
  out$cycles = lapply(which(is.cyc),
                      function(i)
                      {
                        k = which(K[1:ncol(B), i]!=0)
                        v.all = unique(as.vector(ij[k, , drop = FALSE]))
                        sG = graph_from_edgelist(ij[k, , drop = FALSE])
                        tmp.v = v.all[c(1,length(v.all))]
                        p.fwd = shortest_paths(sG, tmp.v[1], tmp.v[2])
                        p.bwd = shortest_paths(sG, tmp.v[2], tmp.v[1])
                        return(node.ix[unique(unlist(c(p.fwd, p.bwd)))])
                      })
  
  out$paths = lapply(which(!is.cyc),
                     function(i)
                     {
                       k = K[1:ncol(B), i]
                       eix = which(k!=0)
                       v.all = unique(as.vector(ij[eix, , drop = FALSE]))
                       sG = graph_from_edgelist(ij[eix, , drop = FALSE])
                       io = B %*% k
                       v.in = Matrix::which(io<0)[1]
                       v.out = Matrix::which(io>0)[1]
                       return(node.ix[unlist(shortest_paths(sG, v.in, v.out))])
                     })
  
  if (length(out$cycles)>0)
  {
    tmp.cix = cbind(unlist(lapply(1:length(out$cycles), function(x) rep(x, length(out$cycles[[x]])))), unlist(out$cycles))
    out$cycles = out$cycles[!duplicated(as.matrix(Matrix::sparseMatrix(tmp.cix[,1], tmp.cix[,2], x = 1)))]
  }
  
  if (length(out$paths)>0)
  {
    tmp.pix = cbind(unlist(lapply(1:length(out$paths), function(x) rep(x, length(out$paths[[x]])))), unlist(out$paths))
    #out$paths = out$paths[!duplicated(as.matrix(Matrix::sparseMatrix(tmp.pix[,1], tmp.pix[,2], x = 1)))]
  }
  
  #if (ALL & length(blank.vertices)>0)
    #out$paths = c(out$paths, lapply(blank.vertices, identity))
  
  return(out)
}

convex.basis = function(A,
                        interval = 80,
                        chunksize = 100,
                        exclude.basis = NULL,
                        exclude.range = NULL,
                        maxchunks = Inf,
                        verbose = T){
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
          
          # remove duplicated rows in H (with respect to sparsity)
          #H = H[!duplicated(as.matrix(H)>ZERO), , drop = F];
          
          # remove rows in H that have subsets in H (with respect to sparsity) ..
          if ((as.numeric(nrow(H))*as.numeric(nrow(H)))>maxchunks){
            print('Exceeding maximum number of chunks in convex.basis computation')
            stop('Exceeding maximum number of chunks in convex.basis computation')
          }
          keep = which(Matrix::colSums(sparse_subset(abs(H)>ZERO, abs(H)>ZERO, chunksize = chunksize, quiet = !verbose))<=1) # <=1 since every H is its own subset
          H = H[keep, , drop = FALSE]
          
          # remove rows in H that have subsets in K_i2
          if (!is.null(K_i2))
            if (nrow(K_i2)>0){
              if ((as.numeric(nrow(K_i2))*as.numeric(nrow(H)))>maxchunks){
                print('Exceeding maximum number of chunks in convex.basis computation')
                stop('Exceeding maximum number of chunks in convex.basis computation')
              }
              keep = which(Matrix::colSums(sparse_subset(abs(K_i2)>ZERO, abs(H)>ZERO, chunksize = chunksize, quiet = !verbose))==0)
              H = H[keep, , drop = FALSE]
            }
          
          # remove rows in H that have subsets in K_i1
          if (!is.null(K_i1))
            if (nrow(K_i1)>0){
              if ((as.numeric(nrow(K_i1))*as.numeric(nrow(H)))>maxchunks)
              {
                print('Exceeding maximum number of chunks in convex.basis computation')
                stop('Exceeding maximum number of chunks in convex.basis computation')
              }
              keep = which(Matrix::colSums(sparse_subset(abs(K_i1)>ZERO, abs(H)>ZERO, chunksize = chunksize, quiet = !verbose))==0)
              H = H[keep, , drop = FALSE]
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

# Set the seed for reproducibility
set.seed(123)

# Create a 15x15 adjacency matrix with random values (0 or 1)
# for a directed graph (asymmetric matrix)
adj_matrix_directed <- matrix(sample(0:1, 15*15, replace = TRUE), nrow = 15, ncol = 15)

# Print the adjacency matrix
adj_matrix_directed

write.csv(adj_matrix_directed,"Scripts/MDM2_Walks_Python/adj_matrix_test.csv")

# create smaller matrix to test all.paths function

adj_matrix_directed_small <- matrix(sample(0:1, 4*4, replace = TRUE), nrow = 4, ncol = 4)

Paths <- all.paths(adj_matrix_directed_small)

Paths

processFile = function(filepath) {
  LineList <- list()
  con = file(filepath, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    LineList <- append(LineList, line)
  }
  
  close(con)
  return(LineList)
}

PythonLists2RLists <- function(StrList){
  NewList <- list()
  for (i in 1:length(StrList)){
    LineString <- StrList[i]
    LineString <- substring(LineString, 3, nchar(LineString)-2)
    LineList <- as.list(strsplit(LineString, "', '"))
    NewList <- append(NewList, LineList)
  }
  for (x in 1:length(NewList)){
    NewList[x] <- list(as.integer(unlist(NewList[x])))
  }
  return(NewList)
}

A1KU19_Cycles <- processFile("Scripts/MDM2_Walks_Python/A1KU19_Cycles.txt")
A1KU19_Cycles_List <- PythonLists2RLists(A1KU19_Cycles)
A1KU19_Cycles_List

for (i in 1:length(NewList)){
  NewList[i] <- as.integer(unlist(NewList[i]))
}

CaseID="A1KU"
CNxt=19
Ggraph <-  gG(jabba=paste0("/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics/Jabba_Calls_TCGA/",CaseID,"/jabba.simple.rds"))
highcopyX = Ggraph[cn>CNxt]        # subset to high copy only

Cycles_gW <- gW(snode.id=A1KU19_Cycles_List, graph=highcopyX, circular = TRUE)

plot(Cycles_gW$gtrack(name = "circular walks"), Cycles_gW$footprint+2e5, cex.label=0.8)
