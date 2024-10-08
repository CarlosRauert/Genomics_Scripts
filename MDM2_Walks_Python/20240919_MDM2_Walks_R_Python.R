library(gGnome)

# Function that can read a File line by Line
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

# Function that turns a list of Python Lists into a list of R Lists
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

# Edited Walks Function that takes Python-generated Walks and Cycles as Input
gGraph$set("public", "walksFromPython", function(field = NULL, greedy = FALSE, verbose = FALSE) {
  
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