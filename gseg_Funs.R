### main functions

##== testing data ==##
# source("fct_SimilarityGraphs.R")
# library(MASS)
# dat1 <- mvrnorm(50, rep(0, 2), Sigma = diag(x = 1, nrow = 2))
# dat2 <- mvrnorm(50, rep(0.5, 2), diag(x = 1, nrow = 2))
# dat <- data.table(group = rep(c("group1", "group2"), c(50, 50)),
#                   rbind(dat1, dat2))
# colnames(dat)[1] <- "group"
# E <- createSimilarityGraph(dat, "MST")
# statistics <- "all"
# n0 <- 10
# n1 <- 10
# n <- nrow(dat)


# single change-point
gseg1 = function(n, E, statistics=c("all","o","w","g","m"), n0=0.05*n, n1=0.95*n, pval.appr=TRUE, skew.corr=TRUE, pval.perm=FALSE, B=100){
  r1 = list()
  
  # range in which to look for CP
  n0 = ceiling(n0)
  n1 = floor(n1)
  
  # Ebynode[[i]]: list of nodes that are connected to i by an edge.
  Ebynode = vector("list", n)
  
  # recode edge list to Ebynode[[i]]
  for(i in 1:n) Ebynode[[i]]=rep(0,0)
  for(i in 1:nrow(E)){
    Ebynode[[E[i,1]]] = c(Ebynode[[E[i,1]]],E[i,2])
    Ebynode[[E[i,2]]] = c(Ebynode[[E[i,2]]],E[i,1])
  }
  
  n0_us = n0
  n1_us = n1
  
  # limit range to [2; n-2] if it is <2, > n-2
  if(n0<2){
    n0=2
  }
  if(n1>(n-2)){
    n1=n-2
  }
  
  # calculate test statistics for all required tests
  r1$scanZ = gcp1bynode(n,Ebynode,statistics,n0,n1)
  
  if (pval.appr==TRUE){
    
    mypval1 = pval1(n,E,Ebynode,r1$scanZ,statistics, skew.corr,n0,n1)
    r1$pval.appr = mypval1
  }
  if (pval.perm==TRUE){
    mypval2 = permpval1(n,Ebynode,r1$scanZ,statistics,B,n0,n1)
    r1$pval.perm = mypval2
  }
  
  if (length(which(!is.na(match(c("o","ori","original","all"),statistics))))>0){
    if(n0_us<=1){
      cat("  Note: Starting index has been set to n0 = 2 as the original edge-count test statistic is not well-defined for t<2. \n")
    }
    if(n1_us>=n-1){
      cat("  Note: Ending index has been set to n1 =", n-2, " as the original edge-count test statistic is not well-defined for t>",n-2,". \n")
    }
  #   cat("Original edge-count scan statistic: \n")
  #   cat("  Estimated change-point location:", r1$scanZ$ori$tauhat, "\n")
  #   cat("  Test statistic:", r1$scanZ$ori$Zmax, "\n")
  #   if (pval.appr==TRUE){
  #     cat("  Approximated p-value:", r1$pval.appr$ori, "\n")
  #   }
  #   if (pval.perm==TRUE){
  #     cat("  p-value from", B, "permutations:", r1$pval.perm$ori$pval, "\n")
  #   }
  # }
  # if (length(which(!is.na(match(c("w","weighted","all"),statistics))))>0){
  #   cat("Weighted edge-count statistic: \n")
  #   if(n0_us<=1){
  #     cat("  Note: Starting index has been set to n0 = 2 as the weighted edge-count test statistic is not well-defined for t<2. \n")
  #   }
  #   if(n1_us>=n-1){
  #     cat("  Note: Ending index has been set to n1 =", n-2, " as the weighted edge-count test statistic is not well-defined for t>",n-2,". \n")
  #   }
  #   cat("  Estimated change-point location:", r1$scanZ$weighted$tauhat, "\n")
  #   cat("  Test statistic:", r1$scanZ$weighted$Zmax, "\n")
  #   if (pval.appr==TRUE){
  #     cat("  Approximated p-value:", r1$pval.appr$weighted, "\n")
  #   }
  #   if (pval.perm==TRUE){
  #     cat("  p-value from", B, "permutations:", r1$pval.perm$weighted$pval, "\n")
  #   }
  # }
  # if (length(which(!is.na(match(c("g","generalized","all"),statistics))))>0){
  #   cat("Generalized edge-count statistic: \n")
  #   if(n0_us<=1){
  #     cat("  Note: Starting index has been set to n0 = 2 as the generalized edge-count test statistic is not well-defined for t<2. \n")
  #   }
  #   if(n1_us>=n-1){
  #     cat("  Note: Ending index has been set to n1 =", n-2, " as the generalized edge-count test statistic is not well-defined for t>",n-2,". \n")
  #   }
  #   cat("  Estimated change-point location:", r1$scanZ$generalized$tauhat, "\n")
  #   cat("  Test statistic:", r1$scanZ$generalized$Zmax, "\n")
  #   if (pval.appr==TRUE){
  #     cat("  Approximated p-value:", r1$pval.appr$generalized, "\n")
  #   }
  #   if (pval.perm==TRUE){
  #     cat("  p-value from", B, "permutations:", r1$pval.perm$generalized$pval, "\n")
  #   }
  # }
  # if (length(which(!is.na(match(c("m","max","all"),statistics))))>0){
  #   cat("Max-type edge-count statistic: \n")
  #   if(n0_us<=1){
  #     cat("  Note: Starting index has been set to n0 = 2 as the max-type edge-count test statistic is not well-defined for t<2. \n")
  #   }
  #   if(n1_us>=n-1){
  #     cat("  Note: Ending index has been set to n1 =", n-2, " as the max-type edge-count test statistic is not well-defined for t>",n-2,". \n")
  #   }
  #   cat("  Estimated change-point location:", r1$scanZ$max.type$tauhat, "\n")
  #   cat("  Test statistic:", r1$scanZ$max.type$Zmax, "\n")
  #   if (pval.appr==TRUE){
  #     cat("  Approximated p-value:", r1$pval.appr$max.type, "\n")
  #   }
  #   if (pval.perm==TRUE){
  #     cat("  p-value from", B, "permutations:", r1$pval.perm$max.type$pval, "\n")
  #   }
  }
  
  return(r1)
}


# the Nu function
Nu = function(x){
  y = x/2
  (1/y)*(pnorm(y)-0.5)/(y*pnorm(y) + dnorm(y))
}

#' gcp1bynode: for single CP
#'   The nodes are numbered by their order in the sequence.
#'   To estimate the change-point, we find the maximum of Z(t), the standardized
#'   version of R(t), between n1 and n2. 
#'
#' @param n no. of nodes
#' @param Ebynode list of nodes connected to node i by an edge (length: n)
#' @param statistics which statistic to calc?
#' @param n0 # lower range limit
#' @param n1 # upper range limit
#'
#' @return scanZ: list of test statistics and underlying processes
gcp1bynode = function(n, Ebynode, statistics="all", n0=ceiling(0.05*n), n1=floor(0.95*n)){
  
  # degree of each node
  nodedeg = rep(0,n)
  for(i in 1:n) nodedeg[i] = length(Ebynode[[i]])
  
  # sum(|G_i|^2)
  sumEisq = sum(nodedeg^2)
  nE = sum(nodedeg)/2
  
  # init 
  g = rep(1,n)   # 1_{i>t}
  R = rep(0,n)   # no of edges connecting obs. before and after t
  R1 = rep(nE,n)  # no. of edges connecting obs. prior to t
  R2 = rep(0,n)  # no. of edges connecting obs. after t
  
  # calculate for each edge R_k(t) (p.7, 2nd paper)
  for(t in 1:(n-1)){
    
    g[t] = 0  # update g (anything <= t is coded as 0, everything after as 1)
    links = Ebynode[[t]]  # all nodes connected to node i
    
    ## first calculate summands for R and R_1
    if(t==1){
      if(length(links)>0){
        
        # sum all linked nodes that are after t (g(t) == 0!)
        R[t] = sum(rep(g[t], length(links)) != g[links])
      } else {
        R[t] = 0
      }
      R1[t]=0  # since no edge possible with both obs. prior to t
      R2[t]=nE-length(links)  # all other edges are after t
      
    } else {  # t > 1
      
      if(length(links)>0){
        # all nodes linked to e that are after e
        add = sum(rep(g[t],length(links)) != g[links])
        # remaining edges that are thus prior to e
        subtract = length(links)-add
        
        # add them to R respectively
        R[t] = R[t-1]+add-subtract
        R1[t] = R1[t-1]+subtract
        
      } else {  # sums remain the same
        R[t] = R[t-1]
        R1[t]=R1[t-1]
      }
    }
    
    ## R2 is all remaining edges
    R2[t] = nE-R[t]-R1[t]
  }
  
  tt = 1:n
  temp=n0:n1
  
  scanZ = list()
  # if we want to calculate original statistic:
  if (length(which(!is.na(match(c("o","ori","original","all"),statistics))))>0){
    # Lemma 2.1, ori paper
    mu.t = nE* 2*tt*(n-tt)/(n*(n-1)) # erwartungswert
    p1.tt = 2*tt*(n-tt)/(n*(n-1))    # p1.tt
    p2.tt = tt*(n-tt)/(n*(n-1))      # 1/2*p1.tt
    p3.tt = 4*tt*(n-tt)*(tt-1)*(n-tt-1)/(n*(n-1)*(n-2)*(n-3))
    # E(X^2) = ...:
    A.tt = p3.tt*nE+(p2.tt-p3.tt)*sumEisq+p3.tt*nE^2
    # calculate Z!
    Z = (mu.t-R)/sqrt(A.tt-mu.t^2)
    Z[n] = 0
    
    tauhat = temp[which.max(Z[n0:n1])]
    ori = list(tauhat = tauhat, Zmax = Z[tauhat], Z = Z, R = R)
    scanZ$ori = ori
  }
  
  # if we want to calculate weighted, max or generalized stat
  if (length(which(!is.na(match(
      c("w","weighted","m","max","g","generalized","all"), statistics))))>0){
    Rw = ((n-tt-1)*R1+(tt-1)*R2)/(n-2)
    mu.Rw = nE*((n-tt-1)*tt*(tt-1)+(tt-1)*(n-tt)*(n-tt-1))/(n*(n-1)*(n-2))
    
    mu.R1 = nE*tt*(tt-1)/(n*(n-1))
    mu.R2 = nE*(n-tt)*(n-tt-1)/(n*(n-1))
    v11 = mu.R1*(1-mu.R1) + 2*(0.5*sumEisq-nE)*(tt*(tt-1)*(tt-2))/(n*(n-1)*(n-2)) + (nE*(nE-1)-2*(0.5*sumEisq-nE))*(tt*(tt-1)*(tt-2)*(tt-3))/(n*(n-1)*(n-2)*(n-3))
    v22 = mu.R2*(1-mu.R2) + 2*(0.5*sumEisq-nE)*((n-tt)*(n-tt-1)*(n-tt-2))/(n*(n-1)*(n-2)) + (nE*(nE-1)-2*(0.5*sumEisq-nE))*((n-tt)*(n-tt-1)*(n-tt-2)*(n-tt-3))/(n*(n-1)*(n-2)*(n-3))
    
    v12 = (nE*(nE-1)-2*(0.5*sumEisq-nE))*tt*(n-tt)*(tt-1)*(n-tt-1)/(n*(n-1)*(n-2)*(n-3)) - mu.R1*mu.R2
    
    var.Rw=((n-tt-1)/(n-2))^2*v11 + 2*((n-tt-1)/(n-2))*((tt-1)/(n-2))*v12+((tt-1)/(n-2))^2*v22
    Zw = -(mu.Rw-Rw)/sqrt(apply(cbind(var.Rw,rep(0,n)),1,max))
    
    if (length(which(!is.na(match(c("w","weighted","m","all"),statistics))))>0){
      tauhat = temp[which.max(Zw[n0:n1])]
      weighted = list(tauhat=tauhat, Zmax=Zw[tauhat], Zw=Zw, Rw=Rw)
      scanZ$weighted = weighted
      
    }
    if (length(which(!is.na(match(c("m","max","g","generalized","all"),statistics))))>0){
      Rd = R1-R2
      Zd = (Rd-(mu.R1-mu.R2))/sqrt(apply(cbind(v11+v22-2*v12,rep(0,n)),1,max))
      
      if (length(which(!is.na(match(c("m","max","all"),statistics))))>0){
        # now makes sure that Inf values aren't counted! 
        M = apply(cbind(abs(Zd), Zw), 1, function(x) suppressWarnings(max(x[is.finite(x)])))
        tauhat = temp[which.max(M[n0:n1])]
        max.type = list(tauhat=tauhat, Zmax=M[tauhat], M=M)
        scanZ$max.type = max.type
      }
      if (length(which(!is.na(match(c("g","generalized","all"),statistics))))>0){
        Z = Zw^2 + Zd^2
        tauhat = temp[which.max(Z[n0:n1])]
        generalized = list(tauhat=tauhat, Zmax=Z[tauhat], S=Z)
        scanZ$generalized = generalized
      }
    }
  }
  return(scanZ)
}


# rho_one = n*h_G
#   p.151 in ori paper! h_1 to h_6 
rho_one = function(n, s, sumE, sumEisq){
  f1 = 4*(n-1)*(2*s*(n-s)-n)
  f2 = ((n+1)*(n-2*s)^2-2*n*(n-1))
  f3 = 4*((n-2*s)^2-n)
  f4 = 4*n*(s-1)*(n-1)*(n-s-1)
  f5 = n*(n-1)*((n-2*s)^2-(n-2))
  f6 = 4*((n-2)*(n-2*s)^2-2*s*(n-s)+n)
  n*(n-1)*(f1*sumE + f2*sumEisq - f3*sumE^2)/(2*s*(n-s)*(f4*sumE + f5*sumEisq - f6*sumE^2))
}

rho_one_Rw = function(n, t){
  -((2*t^2 - 2*n*t + n)*(n^2 - 3*n + 2)^4)/(2*t*(n - 1)^3*(n - 2)^4*(t - 1)*(n^2 - 2*n*t - n + t^2 + t))
}


# p value approximation for single change-point
pval1 = function(n, E, Ebynode, scanZ, statistics="all", skew.corr=TRUE, lower=ceiling(0.05*n), upper=floor(0.95*n)){
  output = list()
  deg = rep(0,n)
  for(i in 1:n) deg[i] = length(Ebynode[[i]])
  sumE = sum(deg)/2
  sumEisq = sum(deg^2)
  
  if (skew.corr==FALSE){
    if (length(which(!is.na(match(c("o","ori","original","all"), statistics))))>0){
      # if(lower<2){
      # lower = 2
      # }
      # if(upper>(n-2)){
      # upper = n-2
      # }
      b = scanZ$ori$Zmax
      if (b>0){
        integrandO = function(s){
          x = rho_one(n,s,sumE,sumEisq)
          x*Nu(sqrt(2*b^2*x))
        }
        pval.ori = dnorm(b)*b*integrate(integrandO, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
      }else{
        pval.ori = 1
      }
      output$ori = min(pval.ori,1)
    }
    
    # if (length(which(!is.na(match(c("w","weighted","m","max","g","generalized","all"),statistics))))>0){
    # if(lower<2){
    # lower = 2
    # }
    # if(upper>(n-2)){
    # upper = n-2
    # }
    # }
    if (length(which(!is.na(match(c("w","weighted","all"), statistics))))>0){
      b = scanZ$weighted$Zmax
      if (b>0){
        integrandW = function(t){
          x = rho_one_Rw(n,t)
          x*Nu(sqrt(2*b^2*x))
        }
        pval.weighted = dnorm(b)*b*integrate(integrandW, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
      }else{
        pval.weighted = 1
      }
      output$weighted = min(pval.weighted,1)
    }
    if (length(which(!is.na(match(c("m","max","all"), statistics))))>0){
      b = scanZ$max.type$Zmax
      if (b>0){
        integrand1 = function(t){
          x1 = n/(2*t*(n - t))
          x1*Nu(sqrt(2*b^2*x1))
        }
        integrand2 = function(t){
          x2 = rho_one_Rw(n,t)
          x2*Nu(sqrt(2*b^2*x2))
        }
        pval_u1 = 2*dnorm(b)*b*integrate(integrand1, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
        pval_u2 = dnorm(b)*b*integrate(integrand2, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
        pval.max.type = as.numeric(1-(1-min(pval_u1,1))*(1-min(pval_u2,1)))
      }else{
        pval.max.type = 1
      }
      output$max.type = pval.max.type
    }
    if (length(which(!is.na(match(c("g","generalized","all"), statistics))))>0){
      b = scanZ$generalized$Zmax
      if (b>0){
        integrandG = function(t,w){
          x1 = n/(2*t*(n - t))
          x2 = rho_one_Rw(n,t)
          2*(x1*cos(w)^2+x2*sin(w)^2)*b*Nu(sqrt(2*b*(x1*cos(w)^2+x2*sin(w)^2)))/(2*pi)
        }
        integrand0 = function(t) {integrate(integrandG,0,2*pi,t=t,subdivisions=3000, stop.on.error=FALSE)$value}
        pval.generalized = dchisq(b,2)*integrate(Vectorize(integrand0),lower,upper,subdivisions=3000, stop.on.error=FALSE)$value
      }else{
        pval.generalized = 1
      }
      output$generalized = min(pval.generalized,1)
    }
    return(output)
  }
  
  x1 = sum(deg*(deg-1))
  x2 = sum(deg*(deg-1)*(deg-2))
  x3 = 0
  for (i in 1:nrow(E)){
    x3 = x3 + (deg[E[i,1]]-1)*(deg[E[i,2]]-1)
  }
  x4 = sum(deg*(deg-1)*(sumE-deg))
  x5 = 0
  for (i in 1:nrow(E)){
    j = E[i,1]
    k = E[i,2]
    x5 = x5 + length(which(!is.na(match(Ebynode[[j]], Ebynode[[k]]))))
  }
  
  if (length(which(!is.na(match(c("o","ori","original","all"), statistics))))>0){
    b = scanZ$ori$Zmax
    if (b>0){
      s = 1:n
      x = rho_one(n,s,sumE,sumEisq)
      p1 = 2*s*(n-s)/(n*(n-1))
      p2 = 4*s*(s-1)*(n-s)*(n-s-1)/(n*(n-1)*(n-2)*(n-3))
      p3 = s*(n-s)*((n-s-1)*(n-s-2) + (s-1)*(s-2))/(n*(n-1)*(n-2)*(n-3))
      p4 = 8*s*(s-1)*(s-2)*(n-s)*(n-s-1)*(n-s-2)/(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5))
      mu = p1*sumE
      sig = sqrt(apply(cbind(p2*sumE + (p1/2-p2)*sumEisq + (p2-p1^2)*sumE^2, rep(0,n)), 1, max))  # sigma
      ER3 = p1*sumE + p1/2*3*x1 + p2*(3*sumE*(sumE-1)-3*x1) + p3*x2 + p2/2*(3*x4-6*x3) + p4*(sumE*(sumE-1)*(sumE-2)-x2-3*x4+6*x3)- 2*p4*x5
      r = (mu^3 + 3*mu*sig^2 - ER3)/sig^3
      theta_b = rep(0,n)
      pos = which(1+2*r*b>0)
      theta_b[pos] = (sqrt((1+2*r*b)[pos])-1)/r[pos]
      ratio = exp((b-theta_b)^2/2 + r*theta_b^3/6)/sqrt(1+r*theta_b)
      a = x*Nu(sqrt(2*b^2*x)) * ratio
      nn = n-length(pos)
      if (nn>0.75*n){
        cat("Extrapolation for skewness-corrected p-value approximation could not be performed. \n")
        cat("Original edge-count statistic: p-value approximation without skewness correction is reported.\n")
        integrand = function(s){
          x = rho_one(n,s,sumE,sumEisq)
          x*Nu(sqrt(2*b^2*x))
        }
        
        pval.ori = dnorm(b)*b*integrate(integrand, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
        output$ori = min(pval.ori,1)
      }else{
        if (nn>=(lower-1)+(n-upper)){
          neg = which(1+2*r*b<=0)
          dif = neg[2:nn]-neg[1:(nn-1)]
          id1 = which.max(dif)
          id2 = id1 + ceiling(0.03*n)
          id3 = id2 + ceiling(0.09*n)
          inc = (a[id3]-a[id2])/(id3-id2)
          a[id2:1] = a[id2+1]-inc*(1:id2)
          a[(n/2+1):n] = a[(n/2):1]
          neg2 = which(a<0)
          a[neg2] = 0
        }
        integrand = function(s){
          a[s]
        }
        
        result = try(dnorm(b)*b*integrate(integrand, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value, silent=T)
        if (is.numeric(result)){
          output$ori = min(result,1)
        }else{
          cat("Original edge-count statistic: p-value approximation without skewness correction is reported.\n")
          b = scanZ$ori$Zmax
          integrand = function(s){
            x = rho_one(n,s,sumE,sumEisq)
            x*Nu(sqrt(2*b^2*x))
          }
          pval.ori = dnorm(b)*b*integrate(integrand, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
          output$ori = min(pval.ori,1)
        }
      }
    }else{
      output$ori = 1
    }
  }
  
  if (length(which(!is.na(match(c("w","weighted","m","max","all"), statistics))))>0){
    t = 1:(n-1)
    A1 = sumE*t*(t-1)/(n*(n-1)) + 3*x1*t*(t-1)*(t-2)/(n*(n-1)*(n-2)) + (3*sumE*(sumE-1)-3*x1)*t*(t-1)*(t-2)*(t-3)/(n*(n-1)*(n-2)*(n-3))  + x2*t*(t-1)*(t-2)*(t-3)/(n*(n-1)*(n-2)*(n-3)) + (6*x3 - 6*x5)*(t*(t-1)*(t-2)*(t-3))/(n*(n-1)*(n-2)*(n-3)) +
      2*x5*(t*(t-1)*(t-2))/(n*(n-1)*(n-2)) + (3*x4+6*x5-12*x3)*t*(t-1)*(t-2)*(t-3)*(t-4)/(n*(n-1)*(n-2)*(n-3)*(n-4)) +
      (sumE*(sumE-1)*(sumE-2)+6*x3-2*x5-x2-3*x4)*t*(t-1)*(t-2)*(t-3)*(t-4)*(t-5)/(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5))
    
    B1 = (sumE*(sumE-1)-x1)*(t*(t-1)*(n-t)*(n-t-1))/(n*(n-1)*(n-2)*(n-3)) + (x4+2*x5-4*x3)*(t*(t-1)*(t-2)*(n-t)*(n-t-1))/(n*(n-1)*(n-2)*(n-3)*(n-4)) +
      (sumE*(sumE-1)*(sumE-2)+6*x3-2*x5-x2-3*x4)*t*(t-1)*(t-2)*(t-3)*(n-t)*(n-t-1)/(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5))
    
    C1 = (sumE*(sumE-1)-x1)*(n-t)*(n-t-1)*t*(t-1)/(n*(n-1)*(n-2)*(n-3)) + (x4+2*x5-4*x3)*(n-t)*(n-t-1)*(n-t-2)*t*(t-1)/(n*(n-1)*(n-2)*(n-3)*(n-4)) +
      (sumE*(sumE-1)*(sumE-2)+6*x3-2*x5-x2-3*x4)*t*(t-1)*(n-t)*(n-t-1)*(n-t-2)*(n-t-3)/(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5))
    
    D1 =  sumE*(n-t)*(n-t-1)/(n*(n-1)) + 3*x1*(n-t)*(n-t-1)*(n-t-2)/(n*(n-1)*(n-2)) + (3*sumE*(sumE-1)-3*x1)*(n-t)*(n-t-1)*(n-t-2)*(n-t-3)/(n*(n-1)*(n-2)*(n-3))  + x2*(n-t)*(n-t-1)*(n-t-2)*(n-t-3)/(n*(n-1)*(n-2)*(n-3)) + (6*x3 - 6*x5)*((n-t)*(n-t-1)*(n-t-2)*(n-t-3))/(n*(n-1)*(n-2)*(n-3)) +
      2*x5*((n-t)*(n-t-1)*(n-t-2))/(n*(n-1)*(n-2)) + (3*x4+6*x5-12*x3)*(n-t)*(n-t-1)*(n-t-2)*(n-t-3)*(n-t-4)/(n*(n-1)*(n-2)*(n-3)*(n-4)) +
      (sumE*(sumE-1)*(sumE-2)+6*x3-2*x5-x2-3*x4)*(n-t)*(n-t-1)*(n-t-2)*(n-t-3)*(n-t-4)*(n-t-5)/(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5))
    
    r1=sumE*(t*(t-1)/(n*(n-1))) + 2*(0.5*sumEisq-sumE)*t*(t-1)*(t-2)/(n*(n-1)*(n-2)) + (sumE*(sumE-1)-(2*(0.5*sumEisq-sumE)))*t*(t-1)*(t-2)*(t-3)/(n*(n-1)*(n-2)*(n-3))
    
    r2=sumE*((n-t)*(n-t-1)/(n*(n-1))) + 2*(0.5*sumEisq-sumE)*(n-t)*(n-t-1)*(n-t-2)/(n*(n-1)*(n-2)) + (sumE*(sumE-1)-(2*(0.5*sumEisq-sumE)))*(n-t)*(n-t-1)*(n-t-2)*(n-t-3)/(n*(n-1)*(n-2)*(n-3))
    
    r12= (sumE*(sumE-1)-(2*(0.5*sumEisq-sumE)))*t*(t-1)*(n-t)*(n-t-1)/(n*(n-1)*(n-2)*(n-3))
    
    t = 1:(n-1)
    x = rho_one_Rw(n,t)
    # for(i in 1:length(x)){
    #   if ((abs(x[i]))=="Inf"){
    #     x[i]=0
    #   }
    # }
    q=(n-t-1)/(n-2)
    p=(t-1)/(n-2)
    
    mu = sumE*(q*t*(t-1)+p*(n-t)*(n-t-1))/(n*(n-1))
    sig1= q^2*r1 + 2*q*p*r12 + p^2*r2 - mu^2
    sig = sqrt(sig1)  # sigma
    ER3 = q^3*A1 + 3*q^2*p*B1 + 3*q*p^2*C1 + p^3*D1
    r =  (ER3- 3*mu*sig^2 - mu^3)/sig^3
    
    b = scanZ$weighted$Zmax
    result.u2 = pval1_sub_2(n,b,r,x,lower,upper)
    r.Rw = r
    x.Rw = x
    
    if (length(which(!is.na(match(c("w","weighted","all"), statistics))))>0){
      if (is.numeric(result.u2) && result.u2 > 0){
        output$weighted = min(result.u2,1)
      }else{
        if (result.u2 ==0){
          cat("Extrapolation for skewness-corrected p-value approximation could not be performed. \n")
        }
        cat("Weighted edge-count statistic: p-value approximation without skewness correction is reported.\n")
        b = scanZ$weighted$Zmax
        if (b>0){
          integrandW = function(t){
            x = rho_one_Rw(n,t)
            x*Nu(sqrt(2*b^2*x))
          }
          pval.weighted = dnorm(b)*b*integrate(integrandW, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
        }else{
          pval.weighted = 1
        }
        output$weighted = min(pval.weighted,1)
      }
    }
    
    if (length(which(!is.na(match(c("m","max","all"), statistics))))>0){
      b = scanZ$max.type$Zmax
      
      t = 1:(n-1)
      x = n/(2*t*(n - t))
      
      q=1
      p=-1
      mu = sumE*(q*t*(t-1)+p*(n-t)*(n-t-1))/(n*(n-1))
      sig1= q^2*r1 + 2*q*p*r12 + p^2*r2 - mu^2
      sig = sqrt(apply(cbind(sig1, rep(0,n-1)), 1, max))  # sigma
      ER3 = q^3*A1 + 3*q^2*p*B1 + 3*q*p^2*C1 + p^3*D1
      r =  (ER3- 3*mu*sig^2 - mu^3)/sig^3
      
      
      result.u1 = pval1_sub_1(n,b,r,x,lower,upper)
      result.u2 = pval1_sub_2(n,b,r.Rw,x.Rw,lower,upper)
      
      if (!is.numeric(result.u1) || !is.numeric(result.u2) || result.u1 ==0 || result.u2 ==0){
        if(result.u1 ==0 || result.u2 == 0){
          cat("Extrapolation for skewness-corrected p-value approximation could not be performed. \n")
        }
        cat("Max-type edge-count statistic: p-value approximation without skewness correction is reported.\n")
        b = scanZ$max.type$Zmax
        if (b>0){
          integrand1 = function(t){
            x1 = n/(2*t*(n - t))
            x1*Nu(sqrt(2*b^2*x1))
          }
          integrand2 = function(t){
            x2 = rho_one_Rw(n,t)
            x2*Nu(sqrt(2*b^2*x2))
          }
          pval_u1 = 2*dnorm(b)*b*integrate(integrand1, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
          pval_u2 = dnorm(b)*b*integrate(integrand2, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
          pval.max.type = as.numeric(1-(1-min(pval_u1,1))*(1-min(pval_u2,1)))
        }else{
          pval.max.type = 1
        }
        output$max.type = pval.max.type
      }else{
        output$max.type = 1-(1-min(result.u1,1))*(1-min(result.u2,1))
      }
    }
  }
  
  # for generalized edge-count test, the approximated p-value without skewness correction is reported
  if (length(which(!is.na(match(c("g","generalized","all"), statistics))))>0){
    b = scanZ$generalized$Zmax
    if (b>0){
      integrandG = function(t,w){
        x1 = n/(2*t*(n - t))
        x2 = rho_one_Rw(n,t)
        2*(x1*cos(w)^2+x2*sin(w)^2)*b*Nu(sqrt(2*b*(x1*cos(w)^2+x2*sin(w)^2)))/(2*pi)
      }
      integrand0 = function(t) {integrate(integrandG,0,2*pi,t=t,subdivisions=3000, stop.on.error=FALSE)$value}
      pval.generalized = dchisq(b,2)*integrate(Vectorize(integrand0),lower,upper,subdivisions=3000, stop.on.error=FALSE)$value
    }else{
      pval.generalized = 1
    }
    output$generalized = min(pval.generalized,1)
  }
  
  
  #if (length(!is.na(match(c("g","generalized","all"), statistics)))>0){
  # cat("Generalized edge-count statistic: p-value approximation without skewness correction is reported.\n")
  #}
  
  return(output)
}

# p-value approximation for single change-point, sub functions
pval1_sub_1 = function(n,b,r,x,lower,upper){
  if (b<0){
    return(1)
  }
  theta_b = rep(0,n-1)
  pos = which(1+2*r*b>0)
  theta_b[pos] = (sqrt((1+2*r*b)[pos])-1)/r[pos]
  for(i in 1:length(theta_b[pos])){
    if (is.na(theta_b[pos][i])==TRUE){
      theta_b[pos][i]=0
    }
  }
  ratio = exp((b-theta_b)^2/2 + r*theta_b^3/6)/sqrt(1+r*theta_b)
  a = x*Nu(sqrt(2*b^2*x)) * ratio
  
  nn.l = ceiling(n/2)-length(which(1+2*r[1:ceiling(n/2)]*b>0))
  nn.r = ceiling(n/2)-length(which(1+2*r[ceiling(n/2):(n-1)]*b>0))
  if (nn.l>0.35*n || nn.r>0.35*n){
    return(0)
  }
  if (nn.l>=lower){
    neg = which(1+2*r[1:ceiling(n/2)]*b<=0)
    dif = c(diff(neg),n/2-nn.l)
    id1 = which.max(dif)
    id2 = id1 + ceiling(0.03*n)
    id3 = id2 + ceiling(0.09*n)
    inc = (a[id3]-a[id2])/(id3-id2)
    a[id2:1] = a[id2+1]-inc*(1:id2)
  }
  if (nn.r>=(n-upper)){
    neg = which(1+2*r[ceiling(n/2):(n-1)]*b<=0 )
    id1 = min(neg+ceiling(n/2)-1,ceiling(n/2)-1)
    id2 = id1 - ceiling(0.03*n)
    id3 = id2 - ceiling(0.09*n)
    inc = (ratio[id3]-ratio[id2])/(id3-id2)
    ratio[id2:(n-1)] = ratio[id2-1]+inc*((id2:(n-1))-id2)
    ratio[ratio<0]=0
    a[(n/2):(n-1)] = (x*Nu(sqrt(2*b^2*x)) * ratio)[(n/2):(n-1)] # update a after extrapolation 
  }
  neg2 = which(a<0)
  a[neg2] = 0
  integrand = function(s){
    a[s]
  }
  result = try(2*dnorm(b)*b*integrate(integrand, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value, silent=T)
  return(result)
  
}

pval1_sub_2 = function(n,b,r,x,lower,upper){
  if (b<0){
    return(1)
  }
  theta_b = rep(0,n-1)
  pos = which(1+2*r*b>0)
  theta_b[pos] = (sqrt((1+2*r*b)[pos])-1)/r[pos]
  ratio = exp((b-theta_b)^2/2 + r*theta_b^3/6)/sqrt(1+r*theta_b)
  a = x*Nu(sqrt(2*b^2*x)) * ratio
  a_na = which(is.na(a)==TRUE )
  a[a_na] = 0
  nn = n-1-length(pos)
  if (nn>0.75*n){
    return(0)
  }
  if (nn>=(lower-1)+(n-upper)){
    neg = which(1+2*r*b<=0)
    dif = neg[2:nn]-neg[1:(nn-1)]
    id1 = which.max(dif)
    id2 = id1 + ceiling(0.03*n)
    id3 = id2 + ceiling(0.09*n)
    inc = (a[id3]-a[id2])/(id3-id2)
    a[id2:1] = a[id2+1]-inc*(1:id2)
    a[(n/2+1):n] = a[(n/2):1]
    neg2 = which(a<0 | is.na(a)==TRUE)
    a[neg2] = 0
  }
  integrand = function(s){
    a[s]
  }
  result = try(dnorm(b)*b*integrate(integrand, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value, silent=T)
  return(result)
}


#' permpval1: p value from permutation for single change point
#'
#' @param n          num obs.
#' @param Ebynode    list of nodes connected to node i by an edge (length: n)
#' @param scanZ      output of processes by gcp1bynode      
#' @param statistics which CP statistics should I compute
#' @param B          number of permutations
#' @param n0         start of CP region
#' @param n1         end of CP region
#'
#' @return
permpval1 = function(n, Ebynode, scanZ, statistics="all", B=100, n0=ceiling(0.05*n), n1=floor(0.95*n)){
  # Computes the pvalue P(max_{n1<=t<=n2} Z(t)>b) by permuting the nodes in the graph.
  Z.ori = Z.weighted = Z.max.type = Z.generalized = matrix(0,B,n)
  
  ### run through all permutations
  for(b in 1:B){
    if(b%%1000 ==0) {
      cat(b, "permutations completed.\n")
    }
    # permutation of the time indices
    perm = sample(n)
    
    permmatch = rep(0,n)
    for(i in 1:n) permmatch[perm[i]] = i
    
    Ebnstar =  vector("list", n)
    for(i in 1:n){
      oldlinks = Ebynode[[permmatch[i]]]
      Ebnstar[[i]] = perm[oldlinks]  # list of nodes connected to i after permutation
    }
    gcpstar=gcp1bynode(n,Ebnstar,statistics,n0,n1)  # statistics for permuted graph
    
    if (length(which(!is.na(match(c("o","ori","original","all"), statistics))))>0){
      Z.ori[b,] = gcpstar$ori$Z
    }
    if (length(which(!is.na(match(c("w","weighted","all"), statistics))))>0){
      Z.weighted[b,] = gcpstar$weighted$Zw
    }
    if (length(which(!is.na(match(c("m","max","all"), statistics))))>0){
      Z.max.type[b,] = gcpstar$max.type$M
    }
    if (length(which(!is.na(match(c("g","generalized","all"), statistics))))>0){
      Z.generalized[b,] = gcpstar$generalized$Z
    }
  }
  ####
  
  output = list()
  p=1-(0:(B-1))/B  # vector of length n seq from 1 to (1 - B-1/B)
                   # just for plotting purposes I think
  
  if (length(which(!is.na(match(c("o","ori","original","all"), statistics))))>0){
    # if((n0<=1 & n1>=(n-2)) | (n0<=2 & n1>=(n-1))){
    # n0 = 2
    # n1 = n-2
    # }
    if(is.matrix(Z.ori[,n0:n1])){
      maxZ = apply(Z.ori[,n0:n1],1,max)
    } else { # if you test only one CP
      maxZ <- Z.ori[,n0:n1]
    }
    maxZs = sort(maxZ)
    output$ori = list(pval=length(which(maxZs>=scanZ$ori$Zmax))/B,  # how many test statistics were more extreme than the original one?
                      curve=cbind(maxZs,p),
                      maxZs=maxZs, Z=Z.ori)
  }
  
  
  # if (length(which(!is.na(match(c("w","weighted","m","max","g","generalized","all"),statistics))))>0){
  # if((n0<=1 & n1>=(n-2)) | (n0<=2 & n1>=(n-1))){
  # n0 = 2
  # n1 = n-2
  # }
  # }
  if (length(which(!is.na(match(c("w","weighted","all"), statistics))))>0){
    
    if(is.matrix(Z.weighted[,n0:n1])){
      maxZ = apply(Z.weighted[,n0:n1],1,max)
    } else { # if you test only one CP
      maxZ <- Z.weighted[,n0:n1]
    }
    maxZs = sort(maxZ)
    output$weighted = list(pval=length(which(maxZs>=scanZ$weighted$Zmax))/B, curve=cbind(maxZs,p), maxZs=maxZs, Z=Z.weighted)
  }
  
  if (length(which(!is.na(match(c("m","max","all"), statistics))))>0){
    if(is.matrix(Z.max.type[,n0:n1])){
      maxZ = apply(Z.max.type[,n0:n1],1,max)
    } else { # if you test only one CP
      maxZ <- Z.max.type[,n0:n1]
    }
    maxZs = sort(maxZ)
    output$max.type = list(pval=length(which(maxZs>=scanZ$max.type$Zmax))/B, curve=cbind(maxZs,p), maxZs=maxZs, Z=Z.max.type)
  }
  
  if (length(which(!is.na(match(c("g","generalized","all"), statistics))))>0){
    if(is.matrix(Z.generalized[,n0:n1])){
      maxZ = apply(Z.generalized[,n0:n1],1,max)
    } else { # if you test only one CP
      maxZ <- Z.generalized[,n0:n1]
    }
    maxZs = sort(maxZ)
    output$generalized = list(pval=length(which(maxZs>=scanZ$generalized$Zmax))/B, curve=cbind(maxZs,p), maxZs=maxZs, Z=Z.generalized)
  }
  
  return(output)
}