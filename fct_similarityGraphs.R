# requires
library(ade4) # ade4::mstree
library(cccd) # cccd::nng
library(nbpMatching) # nbpMatching::nonbimatch / ::distancematrix
library(data.table); library(ggplot2); # library(gSeg)
library(ggplot2)

#' Title
#'
#' @param data 
#' @param type 
#' @param k 
#' @param dist_measure c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowsky")
#'
#' @return
createSimilarityGraph <- function(data = NULL, type = c("MST", "MDP", "NNG"), k = 1,
                                  dist_measure = "euclidean", distmat = NULL){
    if(type == "MST"){
    if(is.null(distmat)){
      dists <- dist(data[,-1])
    } else {
      dists <- as.dist(distmat)
    }
      
    MST <- ade4::mstree(dists, ngmax = k)  # k-MST
    return(MST)
  }
  
  if(type == "MDP"){
    
    if(is.null(distmat)){
      adjmat <- as.matrix(dist(data[,-1], diag = TRUE, upper =TRUE,
                               method = dist_measure))
      dists <- nbpMatching::distancematrix(adjmat)
    } else {
      dists <- nbpMatching::distancematrix(as.matrix(distmat))
    }
    
    MDP <- nbpMatching::nonbimatch(dists)
    return(MDP)
  }
  
  if(type == "NNG"){
    if(is.null(distmat)){
      
      dists <- as.matrix(dist(data[,-1], diag = TRUE, upper =TRUE,
                               method = dist_measure))
    } else {
      dists <- distmat
    }
    NNG <- cccd::nng(dx = dists, k = 1)
    return(NNG)
  }
}


#' Title
#'
#' @param data must include a group column with group1 and group2
#' @param graph 
#' @param pch_group1 
#' @param pch_group2 
#'
#' @return currently pretty ugly ggplot object
plotSimilarityGraph <- function(data, graph,
                                pch_group1 = 1, pch_group2 = 17){
  dt <- copy(data)
  setDT(dt)
  
  if(!any(colnames(dt) == "group")){
    message("data needs a group column with the values group1 and group2!")
    return()
  }

  # MST
  if(class(graph) == "neig"){
    MST_edges <- unclass(graph)
      
    startPts <- dt[,!"group"][MST_edges[,1]]
    endPts   <- dt[,!"group"][MST_edges[,2]]
    
    linePts <- cbind(startPts, endPts)
    colnames(linePts) <- c("x0", "y0", "x1", "y1")
  }
  
  # MDP
  if(class(graph) == "nonbimatch"){
    MDP_edges <- graph$matches[, c(2,4)]
    
    MDP_edges <- data.table(t(apply(MDP_edges, 1, sort)))
    MDP_edges <- as.matrix(unique(MDP_edges))
    
    startPts <- dt[MDP_edges[,1], 2:3]
    endPts   <- dt[MDP_edges[,2], 2:3]
    
    linePts <- cbind(startPts, endPts)
    colnames(linePts) <- c("x0", "y0", "x1", "y1")
  }
  
  # NNG
  if(class(graph) == "igraph"){
    NNG_edges <- igraph::as_edgelist(graph)
    
    NNG_edges <- data.table(t(apply(NNG_edges, 1, sort)))
    NNG_edges <- as.matrix(unique(NNG_edges))
    
    startPts <- dt[NNG_edges[,1], 2:3]
    endPts   <- dt[NNG_edges[,2], 2:3]
    
    linePts <- cbind(startPts, endPts)
    colnames(linePts) <- c("x0", "y0", "x1", "y1")
  }
  
  sim_graph <- ggplot(aes(x = base::get(colnames(dt)[2]),
                          y = base::get(colnames(dt)[3])), data = dt) +
    geom_point(aes(shape = group), size = 3) +
    geom_segment(data = linePts, aes(x = x0, y = y0, xend = x1, yend = y1)) +
    theme_bw()
  
  return(sim_graph)
}


runGsegTest <- function(data = NULL, type, k = 1, statistics = "all",
                        n0=0.05*nrow(data), n1=0.95*nrow(data),
                        pval.appr = TRUE, pval.perm = !pval.appr,
                        B = 100, distmat = NULL){
  # create similarity graph
  sim_graph <- createSimilarityGraph(data, type, k, distmat = distmat)
  if(is.null(data)) n <- dim(distmat)[1] else n <- nrow(data)
  
  # depending on graph type, read out edges and then run CP test
  if(type == "MST"){
    test_out <- gseg1(n = n, 
                      E = sim_graph,
                      statistics = statistics,
                      n0 = n0, n1 = n1,
                      pval.appr = pval.appr,
                      pval.perm = pval.perm, B = B)
  }
  
  if(type == "MDP"){

    MDP_edges <- sim_graph$matches[, c(2,4)]
    MDP_edges <- data.table(t(apply(MDP_edges, 1, sort)))
    MDP_edges <- as.matrix(unique(MDP_edges))
    E <- ade4::neig(edges = MDP_edges)
    
    # MDP graphs can only compute the weighted and original statistic!
    # Z_w is not defined for them because the variance of R_w is zero!! 
    test_out <- gseg1(n = nrow(data),
                      E = E,
                      statistics = statistics,
                      n0 = n0, n1 = n1,
                      pval.appr = pval.appr,
                      pval.perm = pval.perm, B = B)
  }
  
  if(type == "NNG"){
    NNG_edges <- igraph::as_edgelist(sim_graph)
    NNG_edges <- data.table(t(apply(NNG_edges, 1, sort)))
    NNG_edges <- as.matrix(unique(NNG_edges))
    E <- ade4::neig(edges = NNG_edges)
    
    test_out <- gseg1(n = nrow(data),
                      E = E,
                      statistics = statistics,
                      n0 = n0, n1 = n1,
                      pval.appr = pval.appr,
                      pval.perm = pval.perm, B = B)
  }
  
  return(test_out)
}
