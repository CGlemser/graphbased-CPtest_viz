library(MASS)
library(data.table)
library(compositions)  # for mv lognormal
library(ade4) # ade4::mstree
library(cccd) # cccd::nng
library(nbpMatching) # nbpMatching::nonbimatch / ::distancematrix
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

createDataSet <- function(n, delta, delta_var, tau, d = 2,
                          distribution = c("normal", "lognormal")){
  if(is.null(tau)){
    tau <- floor(0.5*n)
  }
  
  # how many data points are in group1 and group2
  n_group1 <- tau
  n_group2 <- n - n_group1
  Sigma <- diag(x = 1, nrow = d)
  mu <- rep(0, d)
  
  dist <- match.arg(distribution)
  if(dist == "normal"){
    dat1 <- mvrnorm(n_group1, mu, Sigma)
    dat2 <- mvrnorm(n_group2, mu + delta/sqrt(d), Sigma + delta_var/sqrt(d))
    
    dat <- data.table(group = rep(c("group1", "group2"), c(n_group1, n_group2)),
                      rbind(dat1, dat2))
  }
  if(dist == "lognormal"){
    # dat2 <- MASS::mvrnorm(n_group2, mu + delta, Sigma)
    dat1 <- compositions::rlnorm.rplus(n_group1, mu, Sigma)
    dat2 <- compositions::rlnorm.rplus(n_group2, mu + delta/sqrt(d), Sigma)
    
    dat <- data.table(group = rep(c("group1", "group2"), c(n_group1, n_group2)),
                      rbind(dat1, dat2))
  }
  return(dat)
}

extractEdges <- function(graph){
  if(class(graph) == "neig"){
    edges <- unclass(graph)
  }
  
  if(class(graph) == "nonbimatch"){
    edges <- graph$matches[, c(2,4)]
    
    edges <- data.table(t(apply(edges, 1, sort)))
    edges <- as.matrix(unique(edges))
  }
  
  if(class(graph) == "igraph"){
    edges <- igraph::as_edgelist(graph)
    
    edges <- data.table(t(apply(edges, 1, sort)))
    edges <- as.matrix(unique(edges))
  }
  # sort by first, then second column
  edges <- edges[order(edges[,1], edges[,2]),]
  return(edges)
}

# boolean matrix: which edges contribute to R0, R1, R2
find_Redges <- function(E,t){
  E_bool <- (E > t)
  E_bool <- cbind(E_bool, E_bool[,1] != E_bool[,2])  # J = 0, R_0
  E_bool <- cbind(E_bool, E_bool[,1] == E_bool[,2] & E_bool[,1] == FALSE)
  E_bool <- cbind(E_bool, E_bool[,1] == E_bool[,2] & E_bool[,1] == TRUE)
  colnames(E_bool) <- c("gi", "gj", "J0", "J1", "J2")
  
  return(E_bool)
}

#' Title
#'
#' @param data must include a group column with group1 and group2
#' @param graph 
#' @param pch_group1 
#' @param pch_group2 
#' @param edge_highlights boolean vector as long as #edges indicating which ones to print bold 
#'
#' @return currently pretty ugly ggplot object
plotSimilarityGraph_app <- function(graph, dt, edge_highlights = NULL,
                                    vertex_size = 2, edge_size = .5){

  sim_igraph <- igraph::graph_from_edgelist(extractEdges(graph), directed = FALSE)
  
  set.seed(1) # set.seed so layout is always the same!
  layout <- as.data.frame(igraph::layout.fruchterman.reingold(sim_igraph))
  layout$ID <- factor(1:nrow(dt))
  layout$group <- factor(dt$group)
  
  edgelist <- get.data.frame(sim_igraph)
  # order it the same way as extractEdges, so that edge_highlights match!
  edgelist <- edgelist[order(edgelist$from, edgelist$to),]
  edgelist$from.x <- layout$V1[match(edgelist$from, layout$ID)]
  edgelist$from.y <- layout$V2[match(edgelist$from, layout$ID)]
  edgelist$to.x   <- layout$V1[match(edgelist$to, layout$ID)]
  edgelist$to.y   <- layout$V2[match(edgelist$to, layout$ID)]
  
  edgelist$col <- "#2C3E50" 
  
  if(!is.null(edge_highlights)){
    edgelist$lwd <- 0.5*edge_highlights+1
    edgelist$col[edge_highlights] <- "#E74C3C"
  } else {
    edgelist$lwd <- 1
  }
  
  p <- ggplot(data = layout, aes(x = V1, y = V2, col = group)) +
    geom_segment(data = edgelist, aes(x = from.x, xend = to.x, y = from.y, yend = to.y),
                 size = edgelist$lwd, col = edgelist$col) +
    geom_point(aes(shape = group), size = vertex_size + 1, col = "black") +
    geom_point(aes(shape = group), size = vertex_size) +
    scale_color_viridis_d() +
    theme_void(base_size = 20)
  
  return(p)
}


# plotSimilarityGraph_app <- function(data, graph,
#                                 pch_group1 = 1, pch_group2 = 17,
#                                 edge_highlights = NULL){
#   
#   dt <- copy(data)
#   setDT(dt)
#   
#   if(!any(colnames(dt) == "group")){
#     message("data needs a group column with the values group1 and group2!")
#     return()
#   }
#   
#   edges <- extractEdges(graph)
#   
#   startPts <- dt[,!"group"][edges[,1]]
#   endPts   <- dt[,!"group"][edges[,2]]
#     
#   linePts <- cbind(startPts, endPts)
#   colnames(linePts) <- c("x0", "y0", "x1", "y1")
#   linePts$col <- "#2C3E50" 
#   
#   if(!is.null(edge_highlights)){
#     linePts$lwd <- 0.5*edge_highlights+1
#     linePts$col[edge_highlights] <- "#E74C3C"
#   } else {
#     linePts$lwd <- 1
#   }
# 
#   
#   sim_graph <- ggplot(aes(x = base::get(colnames(dt)[2]),
#                           y = base::get(colnames(dt)[3])), data = dt) +
#     geom_point(aes(shape = group, col = group), size = 3) +
#     geom_segment(data = linePts, aes(x = x0, y = y0, xend = x1, yend = y1),
#                  size = linePts$lwd, col = linePts$col) +
#     labs(x = "x", y = "y") +
#     theme_classic(base_size = 20)
#   
#   return(sim_graph)
# }
# 
# n <- 30
# tau <- .5
# 
# dat <- createDataSet(n, 1, tau, "normal")
# graph <- createSimilarityGraph_app(dat, "MDP")
# E <- extractEdges(graph)
# gg_graph <- plotSimilarityGraph(dat, graph)
# 
# E_bool <- find_Redges(E, 15)
# apply(E_bool[,3:5], 2, sum)  # R0, R1, R2

calcExpectations <- function(graph, n, t){
  edges <- extractEdges(graph)
  num_edges <- nrow(edges)
  
  # expectations
  ER0 <- num_edges*2*t*(n-t)/(n*(n-1))
  ER1 <- num_edges*t*(t-1)/(n*(n-1))
  ER2 <- num_edges*(n-t)*(n-t-1)/(n*(n-1))
  
  # sample values
  E_bool <- find_Redges(edges, t)
  R <- apply(E_bool[,3:5], 2, sum)
  
  dat_R <- rbind(c(ER0, ER1, ER2), R)
  dimnames(dat_R) <- list(c("expectation", "sample"), c("R0(t)", "R1(t)", "R2(t)"))
  
  return(dat_R)
}

runGsegTest <- function(data = NULL, sim_graph, statistics = "all",
                        n0=0.05*nrow(data), n1=0.95*nrow(data),
                        pval.appr = FALSE, pval.perm = TRUE, B = 1000){
  # read out edges
  E <- extractEdges(sim_graph)
  # run test
  test_out <- gseg1(n = nrow(data),
                    E = E,
                    statistics = statistics,
                    n0 = n0, n1 = n1,
                    pval.appr = pval.appr,
                    pval.perm = pval.perm, B = B)
  
  # transform test output to table
  mat <- rbind(
    sapply(test_out$scanZ, function(x) x$tauhat),
    sapply(test_out$pval.perm, function(x) x$pval)
  )
  dimnames(mat) <- list(c("tauhat", "p"), c("maxZ", "maxZw", "maxM", "maxS"))
  
  return(mat[,c(1,2,4,3)])
}
