UpdateWs <- function (x, Cs, l1bound, wk = NULL, p = NULL) 
{
  wcss.perfeature <- GetWCSS(x, Cs, wk, p)$wcss.perfeature
  tss.perfeature <- GetWCSS(x, rep(1, nrow(x)))$wcss.perfeature
  lam <- BinarySearch(-wcss.perfeature + tss.perfeature, l1bound)
  ws.unscaled <- soft(-wcss.perfeature + tss.perfeature, lam)
  return(ws.unscaled/l2n(ws.unscaled))
}
BinarySearch <- function (argu, sumabs) 
{
  if (l2n(argu) == 0 || sum(abs(argu/l2n(argu))) <= sumabs) 
    return(0)
  lam1 <- 0
  lam2 <- max(abs(argu)) - 1e-05
  iter <- 1
  while (iter <= 15 && (lam2 - lam1) > (1e-04)) {
    su <- soft(argu, (lam1 + lam2)/2)
    if (sum(abs(su/l2n(su))) < sumabs) {
      lam2 <- (lam1 + lam2)/2
    }
    else {
      lam1 <- (lam1 + lam2)/2
    }
    iter <- iter + 1
  }
  return((lam1 + lam2)/2)
}
l2n <- function (vec) 
{
  return(sqrt(sum(vec^2)))
}
soft <- function (x, d) 
{
  return(sign(x) * pmax(0, abs(x) - d))
}
GetWCSS <- function (x, Cs, wk = NULL, p = NULL, ws = NULL) 
{
  wcss.perfeature <- numeric(ncol(x))
  for (k in unique(Cs)) {
    whichers <- (Cs == k)
    if (sum(whichers) > 1) 
      if (!is.null(wk))
        wcss.perfeature <- wcss.perfeature + apply((scale(x[whichers, 
                                                         ], center = TRUE, scale = FALSE)^2), 2, sum)*(wk[k])^p
      if (is.null(wk))
        wcss.perfeature <- wcss.perfeature + apply(scale(x[whichers, 
                                                           ], center = TRUE, scale = FALSE)^2, 2, sum)
  
  }
  bcss.perfeature <- apply(scale(x, center = TRUE, scale = FALSE)^2, 
                           2, sum) - wcss.perfeature
  if (!is.null(ws)) 
    return(list(wcss.perfeature = wcss.perfeature, wcss = sum(wcss.perfeature), 
                wcss.ws = sum(wcss.perfeature * ws), bcss.perfeature = bcss.perfeature))
  if (is.null(ws)) 
    return(list(wcss.perfeature = wcss.perfeature, wcss = sum(wcss.perfeature), 
                bcss.perfeature = bcss.perfeature))
}