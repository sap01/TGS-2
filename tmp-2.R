#######################################################################################################
## Goal: Learn CLR1.3 net.
## TODO(sap)
##
LearnClr1dot3NetMfi <- function(mut.info.matrix, mi.net.adj.matrix, num.nodes, max.fanin, 
                                output.dirname, init.path) {
  ##------------------------------------------------------------
  ## Begin: Load the Required Libraries
  ##------------------------------------------------------------
  # none
  ##------------------------------------------------------------
  ## End: Load the Required Libraries
  ##------------------------------------------------------------
  
  ## Load external functions
  source(paste(init.path, 'learn_clr.R', sep = '/'))
  
  # dyn.load(paste(init.path, 'src/clr.so', sep = '/'))
  if (.Platform$OS.type == 'windows') {
    ## TODO(sap)
    # dyn.load(paste(init.path, 'src/windows/clr.dll', sep = '/'))
  } else if (.Platform$OS.type == 'unix') {
    dyn.load(paste(init.path, 'src/unix/clr.so', sep = '/'))
  }
  
  ## source(paste(init.path, 'learn_clr.R', sep = '/'))
  mi.net.adj.matrix.wt <- LearnClr(mut.info.matrix, 'CLR1.3') # weighted adj matrix
  
  # Replace 'NaN' with zero. 'NaN' is produced when a corr. variable has variance zero.
  mi.net.adj.matrix.wt[is.nan(mi.net.adj.matrix.wt)] <- 0
  
  # writeLines('\n mi.net.adj.matrix.wt = \n')
  # print(mi.net.adj.matrix.wt)
  save(mi.net.adj.matrix.wt, file = paste(output.dirname, 'mi.net.adj.matrix.wt.RData', sep = '/'))
  
  # For each target node
  for (col.idx in 1:num.nodes) {
    # Weights of the edges with the target node
    edge.wts <- mi.net.adj.matrix.wt[, col.idx]
    
    # Count number of neighbours having positive edge weight
    num.nbrs <- length(edge.wts[edge.wts > 0])
    
    if (num.nbrs >= max.fanin) {
      # Return indices of the top 'max.fanin' number of neighbours w.r.t. edge weight.
      # Tie is broken in favour of the neighbour having smaller index.
      valid.nbrs <- sort(edge.wts, decreasing = TRUE, index.return = TRUE)$ix[1:max.fanin]
      
      mi.net.adj.matrix[valid.nbrs, col.idx] <- 1
      
      ## The following line is not required since 'mi.net.adj.matrix' is initialized
      ## with all zeroes
      # mi.net.adj.matrix[-(valid.nbrs), col.idx] <- 0
    } else if (num.nbrs < max.fanin) {
      # Retain all the neighbours
      mi.net.adj.matrix[edge.wts > 0, col.idx] <- 1
    }
  }
  
  return(mi.net.adj.matrix)
}
############################################################################################