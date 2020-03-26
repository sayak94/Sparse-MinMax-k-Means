#matequal <- function(x, y)
#is.matrix(x) && is.matrix(y) && dim(x) == dim(y) && all(x == y)


minmaxkmeans <- function(F, lam, p_max=0.5, pstep=0.01, beta=0){
  if(length(lam) == 1L){
    M <- lam
    CENTS <- F[sample.int(as.integer(nrow(F)), M), , drop = FALSE]
    #CENTS <- matrix(rnorm(M*ncol(F),mean=0,sd=1), M, ncol(F))
  }else{
    M <- nrow(lam)
    CENTS <- as.matrix(lam)
  }
  
  
  delta = matrix(data = 0, nrow = nrow(F), ncol = M)
  n = nrow(F)
  p = 0
  D = matrix(data = 0, nrow = n, ncol = M)
  v = vector(mode ="numeric",length = M)
  w = vector(mode ="numeric",length = M)
  wp = vector(mode ="numeric",length = M)
  x = vector(mode ="numeric",length = M)
  cluster = vector(mode ="numeric",length = n)
  flag = 1
  for (i in 1:M){
    w[i] = 1/M
    wp[i] = 1/M
  }
  empty = FALSE
  multi_return <- function() {
    my_list <- list( "centers" = CENTS, "flag" = flag, "cluster" = cluster, "x" = x, "p" = p, "wk" = w)
    return(my_list) 
  }
  niter = 0
  ew = 0
  repeat{
    niter <- niter + 1
   
    for (i in 1:n){
      for (k in 1:M){
        D[i,k] = (w[k]^p)*norm(as.matrix(F[i,] - CENTS[k,]), type = "F")
      }
    }
    U = matrix(data = 0, nrow = nrow(F), ncol = M)
    oldclus <- cluster
    for(i in 1:n){
      nw = which.min(D[i,])
      U[i, nw] = 1
      cluster[i] = nw
    }
    
    for(i in 1:M){
      x[i] = sum(U[,i])
    }
    
    for(i in 1:M){
      if (x[i]<1){
        #change here
        empty = TRUE
        p = p - pstep
        
        U = delta
        w = wp
        break;
      }
    }
    if(p < 0){
      p=0
      flag = 0
      break;
    }
    
    
    for(i in 1:M){
      a = which(U[,i] == 1)
      for(j in 1:ncol(F)){
        CENTS[i,j]=mean(F[a,j])
      }
    }
    if((p<p_max)&&(empty==FALSE)){
      delta = U
      wp = w
      p = p + pstep
    }
    z = 0
    for(k in 1:M){
      v[k]=0
      for(i in 1:n){
        v[k] = v[k] + U[i,k]*norm(as.matrix(F[i,] - CENTS[k,]), type = "F")
      }
      z = z + v[k]^(1/(1-p))
    }
    for(k in 1:M){
      w[k] = beta*w[k] + (1-beta)*v[k]^(1/(1-p))/(z)
    }
    
    prev.ew = ew
    ew = 0
    for (i in 1:n){
      for (k in 1:M){
        D[i,k] = (w[k]^p)*norm(as.matrix(F[i,] - CENTS[k,]), type = "F")
        ew <- ew + D[i,k]
      }
    }
    
    if(niter>1)
      #if(identical(oldclus, cluster)){
      if(abs(ew - prev.ew) <= 1e-04 || niter>=10){
        break
      }
  }
  multi_return()
}