n = 30
cls = Y1
ce = vector(mode ="numeric",length = n)
wb = vector(mode ="numeric",length = n)
for(i in 1:n){
  ce[i] <- classError(km[[i]]$Cs, cls)$errorRate
  wb[i] <- km[[i]]$wbound
}

n=20
ce = vector(mode ="numeric",length = n)
for(i in 1:n){
  km <- SparseCl(F, K=2, wbounds = 1.36)
  ce[i] <- classError(km[[1]]$Cs, class)$errorRate
}

start_time <- Sys.time()
km <- SparseCl(F, K=2, wbounds = 1.62)
end_time <- Sys.time()

end_time - start_time

n=20
ce = vector(mode ="numeric",length = n)
for(i in 1:n){
  cl <- maddkmeans(F, 2)
  ce[i] <- classError(cl$cluster,class)$errorRate
  #km <- KMeansSparseCluster(brain.x, K=5, wbounds = km.perm$bestw)
  #ce[i] <- classError(km[[1]]$Cs, brain.y)$errorRate
}