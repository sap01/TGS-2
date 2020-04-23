#' Calculate correctness
#' 
#' Calculate correctness of the reconstructed
#' network.
NULL

#' Calculate correctness of Predicted Network Adjacency List
#' 
#' Calculate correctness of a predicted directed network 
#' against a given true directed network.
#' 
#' @param pred_net_adj_list Adjaceny list of the predicted
#' directed network. It is a data frame as returned by
#' 'GENIE3::getLinkList()'. The data frame has three
#' columns: regulatoryGene, targetGene and weight.
#' For example, 
#'   regulatoryGene targetGene   weight
#' 1             G3        G10 1.132096
#' 2             G3        G10 1.132096
#' ...
#' 
#' @param true_net_adj_mx Adjacency matrix of the true
#' directed network. 1 = edge, 0 = no edge.
#' 
#' @param print_tr_pos Logical. TRUE = print true
#' positive edges in console, e.g. 'G1 -> G2'. 
#' FALSE (default) = do not print. 
#' 
#' @return R object 'correctness'. It is a 1-by-11 matrix.
#'   The column names represent the correctness 
#'   metrics: 'TP', 'TN', 'FP', 'FN', 'TPR', 'FPR',
#'   'FDR', 'PPV', 'ACC', 'MCC' and 'F1'. The
#'   row stores their corresponding values.
#' 
CalcCorrectnessNetPredDiListVsTrueDiMx <-function(
  pred_net_adj_list, 
  true_net_adj_mx,
  print_tr_pos = FALSE) {
  
  tr_pos <- 0
  tr_neg <- 0
  fl_pos <- 0
  fl_neg <- 0
  
  ## Begin: Calculate True and False Pos
  num_pred_edges <- nrow(pred_net_adj_list)
  
  for (edge_idx in 1:num_pred_edges) {
    src_node <- 
      base::as.character(
        pred_net_adj_list[edge_idx, ]$regulatoryGene)
    
    tgt_node <- 
      base::as.character(
        pred_net_adj_list[edge_idx, ]$targetGene)
    
    if (true_net_adj_mx[src_node, tgt_node] == 1) {
      tr_pos <- (tr_pos + 1)
      
      if (print_tr_pos == TRUE) {
        base::print(
          base::paste(src_node,
                      '->',
                      tgt_node,
                      sep = ' '))
      }
    } else {
      fl_pos <- (fl_pos + 1)
    }
  }
  base::rm(edge_idx)
  
  base::rm(num_pred_edges)
  ## End: Calculate True and False Pos
  
  ## Begin: Calculate False Neg
  num_true_edges <- 
    base::length(
      base::which(true_net_adj_mx == 1))
  
  fl_neg <- (num_true_edges - tr_pos)
  
  base::rm(num_true_edges)
  ## End: Calculate False Neg
  
  ## Begin: Calculate True Neg
  num_true_non_edges <- 
    base::length(
      base::which(true_net_adj_mx == 0))
  
  tr_neg <- (num_true_non_edges - fl_pos)
  
  base::rm(num_true_non_edges)
  ## End: Calculate True Neg
  
  base::rm(pred_net_adj_list, 
           true_net_adj_mx)
  
  #------------------------------------------------------------
  # Begin: Calculate Correctness Metrics
  #------------------------------------------------------------
  TPR <- (tr_pos / (tr_pos + fl_neg))
  FPR <- (fl_pos / (fl_pos + tr_neg))
  
  ## Calculate FDR
  FDR <- NULL
  if ((fl_pos == 0) & (tr_pos == 0)) {
    FDR <- 0
  } else {
    FDR <- (fl_pos / (fl_pos + tr_pos))
  }
  
  ## Calculate PPV
  PPV <- NULL
  if ((fl_pos == 0) & (tr_pos == 0)) {
    PPV <- 0
  } else {
    PPV <- (tr_pos / (tr_pos + fl_pos))
  }
  
  ACC <- ((tr_pos + tr_neg) / 
            (tr_pos + fl_pos + tr_neg + fl_neg))
  
  ## Calculate F1-score
  F1 <- NULL
  if ((PPV == 0) & (TPR == 0)) {
    
    ## Ref: https://github.com/dice-group/gerbil/wiki/Precision,-Recall-and-F1-measure
    F1 <- 0
    
  } else {
    F1 <- (2 * PPV * TPR / (PPV + TPR))
  }
  
  ## Calculate MCC.
  ## '((tr_pos == 0) & (fl_pos == 0))' => Null graph.
  ## '((tr_neg == 0) & (fl_neg == 0))' => Complete graph.
  MCC <- NULL
  if (((tr_pos == 0) & (fl_pos == 0)) | 
      ((tr_neg == 0) & (fl_neg == 0))) {
    
    ## Ref: https://lettier.github.io/posts/2016-08-05-matthews-correlation-coefficient.html
    MCC <- 0
    
  } else {
    MCC <- (((tr_pos * tr_neg) - (fl_neg * fl_pos)) / 
              sqrt((tr_pos + fl_pos) * 
                     (tr_pos + fl_neg) * 
                     (tr_neg + fl_pos) * 
                     (tr_neg+fl_neg)))
  }
  
  ## Calculate AUC under ROC
  # table <- minet::validate(pred_net, true_net)
  # AUC <- minet::auc.roc(table)
  #------------------------------------------------------------
  # End: Calculate Correctness Metrics
  #------------------------------------------------------------
  
  ## Begin: Create the format for correctness
  correctness <- matrix(0, nrow = 1, ncol = 11)
  colnames(correctness) <-
    list('TP',
         'TN',
         'FP',
         'FN',
         'TPR',
         'FPR',
         'FDR',
         'PPV',
         'ACC',
         'MCC',
         'F1')
  ## End: Create the format for correctness
  
  correctness[1, 1] <- tr_pos
  correctness[1, 2] <- tr_neg
  correctness[1, 3] <- fl_pos
  correctness[1, 4] <- fl_neg
  correctness[1, 5] <- TPR
  correctness[1, 6] <- FPR
  correctness[1, 7] <- FDR
  correctness[1, 8] <- PPV
  correctness[1, 9] <- ACC
  correctness[1, 10] <- MCC
  correctness[1, 11] <- F1
  # correctness[1,8] <- AUC
  
  base::return(correctness)
}



#' Calculate correctness of a Directed Network
#'
#' Calculate correctness of a directed network
#' given the number of true positive and false positive
#' edges.
#'
#' @param num_nodes Numerical. Number of nodes.
#' The expected value is a positive integer.
#' 
#' @param num_true_edges Numerical. 
#' Number of true edges.
#' The expected value is a non-negative integer.
#' 
#' @param tr_pos Numerical. 
#' Number of true positive edges.
#' The expected value is a non-negative integer.
#' 
#' @param fl_pos Numerical. 
#' Number of false positive edges.
#' The expected value is a non-negative integer.
#'
#' @return R object 'correctness'. It is a 1-by-11 matrix.
#'   The column names represent the correctness
#'   metrics: 'TP', 'TN', 'FP', 'FN', 'TPR', 'FPR',
#'   'FDR', 'PPV', 'ACC', 'MCC' and 'F1'. The
#'   row stores their corresponding values.
#'
CalcCorrectnessNetPredDiTpFpVsTrueDi <- function(num_nodes,
                                                 num_true_edges,
                                                 tr_pos, 
                                                 fl_pos) {
  tr_neg <- 0
  fl_neg <- 0
  
  ##-----------------------------------------------------------
  ## Begin: Calculate False Neg
  ##-----------------------------------------------------------
  fl_neg <- (num_true_edges - tr_pos)
  ##-----------------------------------------------------------
  ## End: Calculate False Neg
  ##-----------------------------------------------------------
  
  ##-----------------------------------------------------------
  ## Begin: Calculate True Neg
  ##-----------------------------------------------------------
  num_poss_edges <- (num_nodes^2)
  
  num_true_non_edges <- (num_poss_edges - 
                           num_true_edges)
  
  tr_neg <- (num_true_non_edges - fl_pos)
  
  base::rm(num_nodes)
  base::rm(num_true_edges)
  base::rm(num_poss_edges)
  base::rm(num_true_non_edges)
  ##-----------------------------------------------------------
  ## End: Calculate True Neg
  ##-----------------------------------------------------------
  
  
  ##-----------------------------------------------------------
  ## Begin: Calculate Correctness Metrics
  ##-----------------------------------------------------------
  TPR <- (tr_pos / (tr_pos + fl_neg))
  FPR <- (fl_pos / (fl_pos + tr_neg))
  
  ## Calculate FDR
  FDR <- NULL
  if ((fl_pos == 0) & (tr_pos == 0)) {
    FDR <- 0
  } else {
    FDR <- (fl_pos / (fl_pos + tr_pos))
  }
  
  ## Calculate PPV
  PPV <- NULL
  if ((fl_pos == 0) & (tr_pos == 0)) {
    PPV <- 0
  } else {
    PPV <- (tr_pos / (tr_pos + fl_pos))
  }
  
  ACC <- ((tr_pos + tr_neg) /
            (tr_pos + fl_pos + tr_neg + fl_neg))
  
  ## Calculate F1-score
  F1 <- NULL
  if ((PPV == 0) & (TPR == 0)) {
    ## Ref: https://github.com/dice-group/gerbil/wiki/Precision,-Recall-and-F1-measure
    F1 <- 0
    
  } else {
    F1 <- (2 * PPV * TPR / (PPV + TPR))
  }
  
  ## Calculate MCC.
  ## '((tr_pos == 0) & (fl_pos == 0))' => Null graph.
  ## '((tr_neg == 0) & (fl_neg == 0))' => Complete graph.
  MCC <- NULL
  if (((tr_pos == 0) & (fl_pos == 0)) |
      ((tr_neg == 0) & (fl_neg == 0))) {
    ## Ref: https://lettier.github.io/posts/2016-08-05-matthews-correlation-coefficient.html
    MCC <- 0
    
  } else {
    MCC <- (((tr_pos * tr_neg) - (fl_neg * fl_pos)) /
              sqrt((tr_pos + fl_pos) *
                     (tr_pos + fl_neg) *
                     (tr_neg + fl_pos) *
                     (tr_neg + fl_neg)
              ))
  }
  
  ## Calculate AUC under ROC
  # table <- minet::validate(pred_net, true_net)
  # AUC <- minet::auc.roc(table)
  ##-----------------------------------------------------------
  ## End: Calculate Correctness Metrics
  ##-----------------------------------------------------------
  
  ##-----------------------------------------------------------
  ## Begin: Create the format for correctness
  ##-----------------------------------------------------------
  correctness <- matrix(0, nrow = 1, ncol = 11)
  colnames(correctness) <-
    list('TP',
         'TN',
         'FP',
         'FN',
         'TPR',
         'FPR',
         'FDR',
         'PPV',
         'ACC',
         'MCC',
         'F1')
  ##-----------------------------------------------------------
  ## End: Create the format for correctness
  ##-----------------------------------------------------------
  
  correctness[1, 1] <- tr_pos
  correctness[1, 2] <- tr_neg
  correctness[1, 3] <- fl_pos
  correctness[1, 4] <- fl_neg
  correctness[1, 5] <- TPR
  correctness[1, 6] <- FPR
  correctness[1, 7] <- FDR
  correctness[1, 8] <- PPV
  correctness[1, 9] <- ACC
  correctness[1, 10] <- MCC
  correctness[1, 11] <- F1
  # correctness[1,8] <- AUC
  
  base::return(correctness)
  
}
