#######################################################################################################
## Goal: Learn CLR1.2 net. This is exactly same as CLR net.
## Therefore, this function is equivalent to function LearnClrNetMfi().
## The only difference is that in LearnClrNetMfi(), 
## the weighted CLR net adjacency matrix is estimated using the minet
## package's clr() function, whereas in this function, the estimation
## is implemented without depending on any external package.
##
LearnClr1dot2NetMfi <- function(mut.info.matrix, mi.net.adj.matrix, num.nodes, max.fanin, output.dirname) {
  ##------------------------------------------------------------
  ## Begin: Load the Required Libraries
  ##------------------------------------------------------------
  ## none
  ##------------------------------------------------------------
  ## End: Load the Required Libraries
  ##------------------------------------------------------------
  
  ##------------------------------------------------------------
  ## Begin: Estimate weighted CLR net adjacency matrix
  ##------------------------------------------------------------

  
  ## Initialize weighted CLR net adjacency matrix
  # node.names <- rownames(mut.info.matrix)
  # mi.net.adj.matrix.wt <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
  #                                dimnames = c(list(node.names), 
  #                                             list(node.names)))
  mi.net.adj.matrix.wt <- matrix(0, nrow = num.nodes, ncol = num.nodes)
  
  node.mean.sd <- matrix(0, nrow = num.nodes, ncol = 2)
  # rownames(node.mean.sd) <- node.names
  colnames(node.mean.sd) <- c('clr.mean', 'clr.sd')
  
  ## Calculate sample mean and sample standard deviation of each node.
  ## It is calculated acc. to the logic in function 'clr()' in
  ## R package 'minet' (version 3.6.0).
  ## The aforementioned 'clr()' function uses 'clr.cpp' to perform
  ## the calculation. Here, the same logic is re-implemented in R.
  ## Since the in-built 'mean()' and 'sd()' functions in 
  ## R version 3.3.2 follows the exact same logic, therefore, the
  ## re-implementation is straight-forward.
  for (curr.node.idx in 1:num.nodes) {
    ## arithmetic mean
    node.mean.sd[curr.node.idx, 'clr.mean'] <- mean(mut.info.matrix[curr.node.idx, ])
    
    ## var <- 0
    ## for (each sample) {
    ##  sd <- (mean - sample val)
    ##  var <- var + (sd^2)
    ## }
    ## var <- var / (n -1) ## where n = number of samples
    ## sd <- sqrt(var)
    node.mean.sd[curr.node.idx, 'clr.sd'] <- sd(mut.info.matrix[curr.node.idx, ])
  }
  rm(curr.node.idx)
  
  ## Edge weights of the CLR net are calculated acc. to
  ## 'minet::clr()'
  for (tgt.node.idx in 2:num.nodes) {
    for (candidate.parent.idx in 1:(tgt.node.idx - 1)) {
      
      tmp <- 0
      if (node.mean.sd[tgt.node.idx, 'clr.sd'] != 0) {
        tmp <- ((mut.info.matrix[candidate.parent.idx, tgt.node.idx] - 
                  node.mean.sd[tgt.node.idx, 'clr.mean']) / 
          node.mean.sd[tgt.node.idx, 'clr.sd'])
      }
      z.tgt <- max(0, tmp)
      
      tmp <- 0
      if (node.mean.sd[candidate.parent.idx, 'clr.sd'] != 0) {
        tmp <- ((mut.info.matrix[candidate.parent.idx, tgt.node.idx] - 
                  node.mean.sd[candidate.parent.idx, 'clr.mean']) / 
          node.parent.mean.sd[candidate.parent.idx, 'clr.sd'])
      }
      z.parent <- max(0, tmp)
      
      rm(tmp)
      
      clr.edge.wt <- ((z.tgt)^2 + (z.parent)^2)
      clr.edge.wt <- sqrt(clr.edge.wt)
      rm(z.tgt, z.parent)
      
      ## Weighted CLR net adjacency matrix is a symmetric matrix
      mi.net.adj.matrix.wt[candidate.parent.idx, tgt.node.idx] <- clr.edge.wt
      mi.net.adj.matrix.wt[tgt.node.idx, candidate.parent.idx] <- clr.edge.wt
    }
    rm(candidate.parent.idx)
    
  }
  rm(tgt.node.idx)
  ##------------------------------------------------------------
  ## End: Estimate weighted CLR net adjacency matrix
  ##------------------------------------------------------------
  
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