SparseCl <- function (x, K = NULL, wbounds = NULL, nstart = 20, silent = FALSE, 
          maxiter = 6, centers = NULL) 
{
  if (is.null(K) && is.null(centers)) 
    stop("Must provide either K or centers.")
  if (!is.null(K) && !is.null(centers)) {
    if (nrow(centers) != K) 
      stop("If K and centers both are provided, then nrow(centers) must equal K!!!")
    if (nrow(centers) == K) 
      K <- NULL
  }
  if (!is.null(centers) && ncol(centers) != ncol(x)) 
    stop("If centers is provided, then ncol(centers) must equal ncol(x).")
  if (is.null(wbounds)) 
    wbounds <- seq(1.1, sqrt(ncol(x)), len = 20)
  if (min(wbounds) <= 1) 
    stop("wbounds should be greater than 1")
  wbounds <- c(wbounds)
  out <- list()
  if (!is.null(K)) 
    #Cs <- kmeans(x, centers = K, nstart = nstart)$cluster
    mk <- minmaxkmeans(x, K)
    while(!mk$flag)
      mk <- minmaxkmeans(x, K)
    Cs <- mk$cluster
    #Cs <- kmeans(x, centers = cs1)$cluster
  if (is.null(K)) 
    #Cs <- kmeans(x, centers = centers)$cluster
    mk <- minmaxkmeans(x, centers)
    while(!mk$flag)
      mk <- minmaxkmeans(x, centers)
    Cs <- mk$cluster
    #Cs <- kmeans(x, centers = cs1)$cluster
  for (i in 1:length(wbounds)) {
    if (length(wbounds) > 1 && !silent) 
      cat(i, fill = FALSE)
    ws <- rep(1/sqrt(ncol(x)), ncol(x))
    ws.old <- rnorm(ncol(x))
    store.bcss.ws <- NULL
    niter <- 0
    while ((sum(abs(ws - ws.old))/sum(abs(ws.old))) > 1e-04 && 
           niter < maxiter) {
      if (!silent) 
        cat(niter, fill = FALSE)
      niter <- niter + 1
      ws.old <- ws
      if (!is.null(K)) {
        if (niter > 1) 
          mk <- UpdateCs(x, K, ws, Cs)
          Cs <- mk$cluster
      }
      else {
        if (niter > 1) 
          mk <- UpdateCs(x, nrow(centers), ws, Cs)
          Cs <- mk$cluster
      }
      ws <- UpdateWs(x, Cs, wbounds[i], mk$wk, mk$p)
      store.bcss.ws <- c(store.bcss.ws, sum(GetWCSS(x, 
                                                    Cs, mk$wk, mk$p)$bcss.perfeature * ws))
    }
    out[[i]] <- list(ws = ws, Cs = Cs, wcss = GetWCSS(x, 
                                                      Cs, ws=ws, mk$wk, mk$p), crit = store.bcss.ws, wbound = wbounds[i], wk=mk$wk, p=mk$p)
  }
  if (!silent) 
    cat(fill = TRUE)
  class(out) <- "KMeansSparseCluster"
  return(out)
}