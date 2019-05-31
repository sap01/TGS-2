############################################################################################
## Goal: Learn CLR9 net. For each node, retain top 'max.fanin' number of neighbours w.r.t. edge weight
## and remove rest of the edges. Tie is broken in favour of the neighbour having smaller node index.
## If there are less than that number of edges for a node, then retain all its neighbours.
## CLR9 edge weight = (z.parent^2 + z.tgt^2)^0.5 where
## z.parent = max(0, (MI(candidate.parent, tgt) - mean(MI(candidate.parent, .)) / sd(MI(candidate.parent, .)))) and
## z.tgt = max(0, (MI(candidate.parent, tgt) - mean(MI(., tgt)) / sd(MI(., tgt)))).
## If the tgt belongs to time point t_p, then the candidate parent must belong to time point t_(p + 1).
## The MI(*, *) values come from a 'refined' mutual info matrix, where the rows represent 
## the candidate parents and the cols represent the tgts. The refined matrix is derived
## from the raw mutual information matrix between candidate parents and tgts, where 
## the (i, j)^th elt represents the mutual information value between the i^th candidate
## parent and the j^th tgt. The derivation from the raw to refined matrix is done as follows:
## for the i^th candidate parent and the j^th tgt,
## if there exists another candidate parent k (!= i) such that
## MI(i, j) in the raw matrix is less than MI(k, j) in the raw matrix as well as
## less than MI(i, k), which is the mutual info value between the i^th and k^th parent,
## then MI(i, j) is reduced to zero in the refined matrix.
## Otherwise, MI(i, j) is directly copied from the raw matrix to the refined one.
##
LearnClr9NetMfi <- function(input.data.discr.3D, num.nodes, node.names, num.timepts, 
                            max.fanin, mi.net.adj.matrix.list)
{
  ##------------------------------------------------------------
  ## Begin: Load the Required Libraries
  ##------------------------------------------------------------
  ##
  ##------------------------------------------------------------
  ## End: Load the Required Libraries
  ##------------------------------------------------------------
  
  ## Here, each 'time.pt.idx' represents time interval 
  ## ('time.pt.idx', ('time.pt.idx' + 1))
  for (time.pt.idx in 1:(num.timepts - 1)) {
    
    ## Discretized data corr. to the current time interval
    input.data.discr.3D.curr.ival <- input.data.discr.3D[(time.pt.idx:(time.pt.idx + 1)), , ]
    
    candidate.parent.node.names <- c()
    candidate.tgt.node.names <- c()
    for (curr.node.name in node.names) {
      parent.full.name <- paste(curr.node.name, as.character(time.pt.idx), sep = '_t')
      candidate.parent.node.names <- c(candidate.parent.node.names, parent.full.name)
      rm(parent.full.name)
      
      tgt.full.name <- paste(curr.node.name, as.character(time.pt.idx + 1), sep = '_t')
      candidate.tgt.node.names <- c(candidate.tgt.node.names, tgt.full.name)
      rm(tgt.full.name)
    }
    rm(curr.node.name)
    
    ## Initialize raw mutual information matrix between parents and targets with zeroes
    mut.info.matrix.parent.tgt.raw <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                              dimnames = c(list(candidate.parent.node.names), 
                                           list(candidate.tgt.node.names)))
    
    #################################################################################
    # ## Build mutual information matrix
    # for (col.idx in 1:(num.nodes - 1)) {
    #   for (col.idx.2 in (col.idx + 1):num.nodes) {
    #     
    #     ## compute_cmi.R
    #     ## (dim1 == 1) => time.pt.idx
    #     ## (dim1 == 2) => (time.pt.idx + 1)
    #     mut.info <- ComputeCmiPcaCmi(input.data.discr.3D.curr.ival[1, col.idx, ], 
    #                            input.data.discr.3D.curr.ival[2, col.idx.2, ])
    #     
    #     mut.info.matrix.parent.tgt[col.idx, col.idx.2] <- mut.info
    #     mut.info.matrix.parent.tgt[col.idx.2, col.idx] <- mut.info
    #   }
    #   rm(col.idx.2)
    # }
    # rm(col.idx)
    #################################################################################
    ## Build mutual information matrix between parents and targets
    for (parent.idx in 1:num.nodes) {
      for (tgt.idx in 1:num.nodes) {
        
        ## compute_cmi.R
        ## (dim1 == 1) => time.pt.idx
        ## (dim1 == 2) => (time.pt.idx + 1)
        mut.info <- ComputeCmiPcaCmi(input.data.discr.3D.curr.ival[1, parent.idx, ], 
                                     input.data.discr.3D.curr.ival[2, tgt.idx, ])
        
        ## Mutual info matrix is asymmetric in this case
        mut.info.matrix.parent.tgt.raw[parent.idx, tgt.idx] <- mut.info
      }
      rm(tgt.idx)
    }
    rm(parent.idx)
    #################################################################################

    #################################################################################
    ## Build mutual information matrix between parents and parents
    
    ## Initialize mutual information matrix between parents and parents with zeroes
    mut.info.matrix.parent.parent <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                                       dimnames = c(list(candidate.parent.node.names), 
                                                    list(candidate.parent.node.names)))
    
    for (parent.idx in 1:(num.nodes - 1)) {
      for (parent.idx.2 in (parent.idx + 1):num.nodes) {
        
        ## compute_cmi.R
        ## (dim1 == 1) => time.pt.idx
        ## (dim1 == 2) => (time.pt.idx + 1)
        mut.info <- ComputeCmiPcaCmi(input.data.discr.3D.curr.ival[1, parent.idx, ], 
                                     input.data.discr.3D.curr.ival[1, parent.idx.2, ])
        
        ## Mutual info matrix is symmetric in this case
        mut.info.matrix.parent.parent[parent.idx, parent.idx.2] <- mut.info
        mut.info.matrix.parent.parent[parent.idx.2, parent.idx] <- mut.info
      }
      rm(parent.idx.2)
    }
    rm(parent.idx)
    #################################################################################    
    
    #################################################################################
    ## Refine the raw parent-tgt mutual information matrix with the help of 
    ## parent-parent mutual info matrix
    
    ## Initialize the refined matrix with the raw one
    mut.info.matrix.parent.tgt.refined <- mut.info.matrix.parent.tgt.raw
    
    for (parent.idx in 1:num.nodes) {
      for (tgt.idx in 1:num.nodes) {
        parent.tgt.mut.info.raw <- mut.info.matrix.parent.tgt.raw[parent.idx, tgt.idx]
        for (parent.idx.2 in num.nodes) {
          if (parent.idx.2 != parent.idx) {
            parent.tgt.mut.info.raw.2 <- mut.info.matrix.parent.tgt.raw[parent.idx.2, tgt.idx]
            
            parent.parent.mut.info <- mut.info.matrix.parent.parent[parent.idx, parent.idx.2]
            
            if (parent.tgt.mut.info.raw < min(parent.tgt.mut.info.raw.2, parent.parent.mut.info)) {
              mut.info.matrix.parent.tgt.refined[parent.idx, tgt.idx] <- 0
            }
          }
        }
      }
      rm(tgt.idx)
    }
    rm(parent.idx)
    rm(mut.info.matrix.parent.parent, mut.info.matrix.parent.tgt.raw)
    
    #################################################################################
    
        
    ## Initialize unweighted adjacency matrix of the CLR net
    ## corr. to the current time interval
    mi.net.adj.matrix.wt <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                                   dimnames = c(list(candidate.parent.node.names), 
                                                list(candidate.tgt.node.names)))
    
    candidate.parent.mean.sd <- matrix(0, nrow = num.nodes, ncol = 2)
    rownames(candidate.parent.mean.sd) <- candidate.parent.node.names
    colnames(candidate.parent.mean.sd) <- c('clr.mean', 'clr.sd')
    
    ## Calculate sample mean and sample standard deviation of the given nodes.
    ## It is calculated acc. to the logic in function 'clr()' in
    ## R package 'minet' (version 3.36.0).
    ## The aforementioned 'clr()' function uses 'clr.cpp' to perform
    ## the calculation. Here, the same logic is re-implemented in R.
    ## Since the in-built 'mean()' and 'sd()' functions in 
    ## R version 3.3.2 follows the exact same logic, therefore, the
    ## re-implementation is straight-forward.
    for (parent.name in candidate.parent.node.names) {
      ## arithmetic mean
      candidate.parent.mean.sd[parent.name, 'clr.mean'] <- mean(mut.info.matrix.parent.tgt.refined[parent.name, ])
      
      ## var <- 0
      ## for (each sample) {
      ##  sd <- (mean - sample val)
      ##  var <- var + (sd^2)
      ## }
      ## var <- var / (n -1) ## where n = number of samples
      ## sd <- sqrt(var)
      candidate.parent.mean.sd[parent.name, 'clr.sd'] <- sd(mut.info.matrix.parent.tgt.refined[parent.name, ])
    }
    rm(parent.name)
    
    for (tgt.node.name in candidate.tgt.node.names) {
      
      tgt.clr.mean <- mean(mut.info.matrix.parent.tgt.refined[, tgt.node.name])
      tgt.clr.sd <- sd(mut.info.matrix.parent.tgt.refined[, tgt.node.name])
      
      ## Edge weights of the CLR net are calculated acc. to
      ## 'minet::clr()'
      for (candidate.parent.name in candidate.parent.node.names) {
        
        tmp <- 0
        if (tgt.clr.sd != 0) {
          tmp <- (mut.info.matrix.parent.tgt.refined[candidate.parent.name, tgt.node.name] - tgt.clr.mean) / tgt.clr.sd
        }
        z.tgt <- max(0, tmp)
        
        tmp <- 0
        if (candidate.parent.mean.sd[candidate.parent.name, 'clr.sd'] != 0) {
          tmp <- (mut.info.matrix.parent.tgt.refined[candidate.parent.name, tgt.node.name] - 
                    candidate.parent.mean.sd[candidate.parent.name, 'clr.mean']) / 
            candidate.parent.mean.sd[candidate.parent.name, 'clr.sd']
        }
        z.parent <- max(0, tmp)
        
        rm(tmp)
        
        clr.edge.wt <- ((z.tgt)^2 + (z.parent)^2)
        clr.edge.wt <- sqrt(clr.edge.wt)
        rm(z.tgt, z.parent)
        
        mi.net.adj.matrix.wt[candidate.parent.name, tgt.node.name] <- clr.edge.wt
      }
      rm(candidate.parent.name)
      
    }
    rm(tgt.node.name)
    
    ##############################################################
    ## Begin:
    ## Estimate the unweighted adjacency matrix 'mi.net.adj.matrix'
    ## using the weighted adjacency matrix 'mi.net.adj.matrix.wt'
    ## and 'max.fanin'
    ##############################################################
    
    ## Initialize unweighted adjacency matrix of the CLR net
    ## corr. to the current time interval
    mi.net.adj.matrix <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                                dimnames = c(list(candidate.parent.node.names), 
                                             list(candidate.tgt.node.names)))
    
    ## For each target node
    for (col.idx in 1:num.nodes) {
      
      ## Weights of the edges with the target node
      edge.wts <- mi.net.adj.matrix.wt[, col.idx]
      
      ## Count number of neighbours having positive edge weight
      num.nbrs <- length(edge.wts[edge.wts > 0])
      
      if (num.nbrs >= max.fanin) {
        
        ## Return indices of the top 'max.fanin' number of neighbours w.r.t. edge weight.
        ## Tie is broken in favour of the neighbour having smaller index.
        valid.nbrs <- sort(edge.wts, decreasing = TRUE, index.return = TRUE)$ix[1:max.fanin]
        
        mi.net.adj.matrix[valid.nbrs, col.idx] <- 1
        
      } else if (num.nbrs < max.fanin) {
        
        # Retain all the neighbours
        mi.net.adj.matrix[edge.wts > 0, col.idx] <- 1
      }
    }
    rm(col.idx)
    
    ##############################################################
    ## End:
    ## Estimate the unweighted adjacency matrix 'mi.net.adj.matrix'
    ## using the weighted adjacency matrix 'mi.net.adj.matrix.wt'
    ## and 'max.fanin'
    ##############################################################
    
    
    mi.net.adj.matrix.list[[time.pt.idx]] <- mi.net.adj.matrix
  }
  rm(time.pt.idx)
  
  return(mi.net.adj.matrix.list)
}

############################################################################################