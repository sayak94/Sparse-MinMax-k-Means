sp.perm <- function (x, K = NULL, nperms = 25, wbounds = NULL, silent = FALSE, 
                     nvals = 10, centers = NULL) 
{
  if (is.null(wbounds)) 
    wbounds <- exp(seq(log(1.2), log(sqrt(ncol(x)) * 0.9), 
                       len = nvals))
  if (min(wbounds) <= 1) 
    stop("Wbounds should be greater than 1, since otherwise only one weight will be nonzero.")
  if (length(wbounds) < 2) 
    stop("Wbounds should be a vector of at least two elements.")
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
  permx <- list()
  nnonzerows <- NULL
  for (i in 1:nperms) {
    permx[[i]] <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
    for (j in 1:ncol(x)) permx[[i]][, j] <- sample(x[, j])
  }
  tots <- NULL
  out <- SparseCl(x, K, wbounds = wbounds, silent = silent, 
                             centers = centers)
  for (i in 1:length(out)) {
    nnonzerows <- c(nnonzerows, sum(out[[i]]$ws != 0))
    bcss <- GetWCSS(x, out[[i]]$Cs, out[[i]]$wk, out[[i]]$p)$bcss.perfeature
    tots <- c(tots, sum(out[[i]]$ws * bcss))
  }
  permtots <- matrix(NA, nrow = length(wbounds), ncol = nperms)
  for (k in 1:nperms) {
    if (!silent) 
      cat("Permutation ", k, "of ", nperms, fill = TRUE)
    perm.out <- SparseCl(permx[[k]], K, wbounds = wbounds, 
                                    silent = silent, centers = centers)
    for (i in 1:length(perm.out)) {
      perm.bcss <- GetWCSS(permx[[k]], perm.out[[i]]$Cs, perm.out[[i]]$wk, perm.out[[i]]$p)$bcss.perfeature
      permtots[i, k] <- sum(perm.out[[i]]$ws * perm.bcss)
    }
  }
  gaps <- (log(tots) - apply(log(permtots), 1, mean))
  out <- list(tots = tots, permtots = permtots, nnonzerows = nnonzerows, 
              gaps = gaps, sdgaps = apply(log(permtots), 1, sd), wbounds = wbounds, 
              bestw = wbounds[which.max(gaps)])
  if (!silent) 
    cat(fill = TRUE)
  class(out) <- "KMeansSparseCluster.permute"
  return(out)
}