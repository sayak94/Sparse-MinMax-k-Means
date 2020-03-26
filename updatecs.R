UpdateCs <- function (x, K, ws, Cs) 
{
  x <- x[, ws != 0]
  z <- sweep(x, 2, sqrt(ws[ws != 0]), "*")
  nrowz <- nrow(z)
  mus <- NULL
  if (!is.null(Cs)) {
    for (k in unique(Cs)) {
      if (sum(Cs == k) > 1) 
        mus <- rbind(mus, apply(z[Cs == k, ], 2, mean))
      if (sum(Cs == k) == 1) 
        mus <- rbind(mus, z[Cs == k, ])
    }
  }
  if (is.null(mus)) {
    #km <- kmeans(z, centers = K, nstart = 10)
    km <- minmaxkmeans(z, K)
    while(!km$flag)
      km <- minmaxkmeans(z, K)
    #km <- kmeans(z, centers = mk$centers)
  }
  else {
    distmat <- as.matrix(dist(rbind(z, mus)))[1:nrowz, (nrowz + 
                                                          1):(nrowz + K)]
    nearest <- apply(distmat, 1, which.min)
    if (length(unique(nearest)) == K) {
      #km <- kmeans(z, centers = mus)
      km <- minmaxkmeans(z, mus)
      while(!km$flag)
        km <- minmaxkmeans(z, mus)
      #km <- kmeans(z, centers = mk$centers)
    }
    else {
      #km <- kmeans(z, centers = K, nstart = 10)
      km <- minmaxkmeans(z, K)
      while(!km$flag)
        km <- minmaxkmeans(z, K)
      #km <- kmeans(z, centers = mk$centers)
    }
  }
  return(km)
}