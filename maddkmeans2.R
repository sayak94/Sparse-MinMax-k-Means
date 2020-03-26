maddkmeans2 <- function(X, k){
  n = as.integer(nrow(X))
  d = as.integer(ncol(X))
  
  expdist<-function(i, j){
    d <- (-1)*abs(X[i,] - X[j,])
    d <- 1 - exp(d)
    d <- sum(d)
    return(d)
  }
  
  #exp distance matrix
  D = matrix(data = 0, nrow = n, ncol = n)
  md = matrix(data = 0, nrow = n, ncol = n)
  
  mydist<-function(i, j){
    d <- abs(D[i,] - D[j,])
    d[i] = 0
    d[j] = 0
    d <- sum(d)
    d <- 1/(n-2)*d
    return(d)
  }
  
  for(i in 1:n){
    for(j in 1:n){
      if(i != j){
        D[i,j]<-expdist(i, j)
      }
    }
  }
  
  D <- 1/d*D
  
  #calculating rho2 as per paper
  for(i in 1:n){
    for(j in 1:n){
      md[i,j] = mydist(i, j)
    }
  }
  
  cluster = vector(mode = "numeric", length = n)
  
  #initial cluster assignments
  for(i in 1:n){
    cluster[i] = sample(1:k, 1)
  }
  
  niter = 0
  repeat{
    niter <- niter + 1
    U = matrix(data = 0, nrow = n, ncol = k)
    oldclus <- cluster
    
    #u th observation
    for(u in 1:n){
      #j th cluster
      for(j in 1:k){
        
        #list of observations in jth cluster
        cls <- which(cluster %in% j)
        
        #each observation z from cluster j
        for(z in cls){
          U[u,j] <- U[u,j] + md[u,z]^2
        }
        
        U[u,j] <- U[u,j]*1/length(cls)
        
      }
      
      nw = which.min(U[u,])
      cluster[u] = nw
    }
    
    if(niter>1)
      if(identical(oldclus, cluster)){
        break;
      }
  }
  
  multi_return <- function() {
    my_list <- list( "iters" = niter, "cluster" = cluster)
    return(my_list) 
  }
  
  multi_return()
  
}