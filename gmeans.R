


library(MASS)

# set random seed
set.seed(1)

# set alpha value
alpha <- 0.0001
CV <- 1.8692



#F=scale(F)
# generalize for n-dimensional F
p <- ncol(F) + 1

# initialize one cluster membership for all F points
F = cbind(F,rep(1, nrow(F)))

# initalize first center
first <- minmaxkmeans(F[,-p], 1)
centers <- first$centers

# L2 norm calculation
norm_vec <- function(x) sqrt(sum(x^2))

# keep track of number of clusters
numclusters = 1

repeat{
  numcenters <- nrow(centers)
  done = 0
  for(j in 1:numcenters)
  {
    cj <- centers[j, ]
    X <- subset(F, F[,p] == j)
    if(nrow(X) <= 5)
    {
      next
    }
    E <- eigen(cov(as.matrix(X[, -p])))
    m <- E$vectors[, 1]*(sqrt(2*E$values[1]/pi))
    newc1 <- centers[j, ] + m
    newc2 <- centers[j, ] - m
    km1 <- minmaxkmeans(X[, -p], rbind(newc1, newc2))
    km <- kmeans(X[, -p], km1$centers)
    centersdash <- km$centers
    newclusters = km$cluster
    v = centersdash[1, ] - centersdash[2, ]
    Xdash <- as.matrix(X[, -p]) %*% as.matrix(v) / norm_vec(v)
    Xdash = (Xdash - mean(Xdash)) / sd(Xdash)
    sorted <- sort(Xdash)
    z <- pnorm(sorted)
    n = length(z)
    sum = 0
    for(i in 1:n)
    {
      zi = z[i]
      zni = z[n + 1 - i]
      sum = sum + (2 * i - 1) * (log(zi) + log(1 - zni))
    }
    Asq = -n - sum/n
    Asq_st = Asq * (1 + 4 / n - 25 / (n ^ 2))
    if((Asq_st < -CV) | (Asq_st > CV))
    {
      # reject the null hypothesis and accept new centers
      done = 1 #to indicate that we are not yet done
      numclusters = numclusters + 1
      centers[j, ] = centersdash[1, ]
      centers <- rbind(centers, centersdash[2, ])
      # got to be doubly careful about ordering of centers and cluster numbers
      d1 = dist(rbind(X[which(newclusters == 1)[1], -p], centersdash[1, ]))
      d2 = dist(rbind(X[which(newclusters == 1)[1], -p], centersdash[2, ]))
      if(d1 < d2)
      {
        # means that cluster "1" has center centersdash[1]
        newclusters[newclusters == 2] = numclusters
        newclusters[newclusters == 1] = j
      }else{
        # means that cluster "1" has center centersdash[2]
        newclusters[newclusters == 1] = numclusters
        newclusters[newclusters == 2] = j
      }
      F[which(F[,p] == j), p] = newclusters
      
      
    }
  }
  if(done == 0)
  {
    break;
  }
}