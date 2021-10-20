GMDF<-function (object, k, k1, a, lambda.a, output.file, thresh = 1e-04, max.iters = 100, 
                nrep = 1, A.init = NULL,H.init = NULL, W.init = NULL, V.init = NULL,
                rand.seed = 1, print.obj = F,int.val = c(0.1,0.1),param = NULL){
  
  X.update <- function(){
    f<-function(i){
      Xi<-sqrt(a[i,1]) * H[[1]][[i]] %*% A[[1]]
      if(p==1){return(Xi)}
      for(j in 2:p){
        Xi<-Xi + sqrt(a[i,j]) * H[[j]][[i]] %*% A[[j]]
      }
      return(Xi)
    }
    X <- lapply(1:N, f)
    return(X)
  }
  
  W.update<-function(idx = 1:length(X)){
    v<-lapply(1:N,function(i1){E[[i1]] - X[[i1]]})
    W<-solveNNLS(rbindlist(Hw),rbindlist(v))
    return(W)
  }
  
  H.update<-function(i){
    Hw1<-t(solveNNLS(t(W),t(E[[i]] - X[[i]])))
    return(Hw1)
  }
  
  Hij.update<-function(i,j1){
    X.ij <- sqrt(a[i,j1]) * H[[j1]][[i]] %*% A[[j1]]
    Hij <- t(solveNNLS(rbind(sqrt(a[i,j1]) * t(A[[j1]]),
                             sqrt(lambda.a[j1]*a[i,j1]) * t(A[[j1]])),
                       rbind(t(E[[i]] - Hw[[i]] %*% W - (X[[i]] - X.ij)),
                             matrix(0, nrow = g, ncol = ns[i]))))
    return(Hij)
  }
  
  A.update<-function(j1){
    idx<-which(a[,j1]>0)
    f1<-function(i) {
      M <- E[[i]] - Hw[[i]] %*% W - (X[[i]] - sqrt(a[i,j1])* H[[j1]][[i]] %*% A[[j1]])
      return(M)
    }
    
    aH1<-lapply(idx, function(i) return(sqrt(a[i,j1])*H[[j1]][[i]]))
    aH2<-lapply(idx, function(i) return(sqrt(lambda.a[j1]*a[i,j1])*H[[j1]][[i]]))
    aH<-c(aH1,aH2)
    aH<-rbindlist(aH)
    zeroM <- matrix(0, nrow = sum(unlist(ns[idx])),ncol = g)
    M <- rbind(rbindlist(lapply(idx, f1)),zeroM)
    A1 <- solveNNLS(aH,M)
    return(A1)
  }
  
  obj.update<-function(){
    obj1 <- sum(sapply(1:N, function(i) {
      norm(E[[i]] - (Hw[[i]] %*% W + X[[i]]), "F")^2
    }))
    f<-function(j){
      obj2 <- sum(sapply(1:N, function(i) {
        lambda.a[j] * a[i,j] * norm(H[[j]][[i]] %*% A[[j]], "F")^2
      }))
      return(obj2)
    }
    obj2 <- sum(sapply(1:p,f))
    obj <- obj1 + obj2
    return(obj)
  }
  
  object <- removeMissingObs(object, slot.use = "scale.data", 
                             use.cols = F)
  E <- object@scale.data
  # E<-ligerex@raw.data
  N <- length(E)
  ns <- sapply(E, nrow)
  p <- ncol(a)
  if (k >= min(ns)) {
    stop(paste0("Select k lower than the number of cells in smallest dataset: ", 
                min(ns)))
  }
  tmp <- gc()
  g <- ncol(E[[1]])
  if (k >= g) {
    stop(paste0("Select k lower than the number of variable genes:", 
                g))
  }
  W_m <- matrix(0, k, g)
  A_m <- lapply(1:p, function(i) {
    matrix(0, k1, g)
  })
  get.H0 <- function(k){v<-lapply(ns, function(n) {matrix(0, n, k)})}
  
  get.H <- function(k,int.val) {lapply(ns, function(n) {matrix(abs(runif(n * k, 0, int.val)), n, k)})}
  Hw_m <- get.H0(k)
  H_m <- lapply(1:p, function(x) get.H0(k1))
  
  tmp <- gc()
  best_obj <- Inf
  run_stats <- matrix(0, nrow = nrep, ncol = 2)
  objAA <- matrix(nrow = max.iters,ncol = nrep)
  for (i in 1:nrep) {
    set.seed(rand.seed + i - 1)
    start_time <- Sys.time()
    W <- matrix(abs(runif(g * k, 0, int.val)), k, g)
    A <- lapply(1:p, function(i) {
      matrix(abs(runif(g * k1, 0, int.val[2])), k1, g)
    })
    Hw <- get.H(k,int.val[1])
    H <- lapply(1:p, function(x) get.H(k1,int.val[2]))
    
    tmp <- gc()
    if (!is.null(W.init)) {W <- W.init}
    if (!is.null(H.init)) {H <- H.init}
    if (!is.null(A.init)) {A <- A.init}
    
    delta <- 1
    iters <- 0
    pb <- txtProgressBar(min = 0, max = max.iters, style = 3)
    X <- X.update()
    obj0 <- obj.update()
    tmp <- gc()
    
    objA<-c(0)
    while (delta > thresh & iters < max.iters) {
      for(i1 in 1:N){
        for(j1 in 1:p){
          H[[j1]][[i1]]<-Hij.update(i1,j1 = j1)
          X <- X.update()
        }
      }
      A <- lapply(1:p, A.update); tmp <- gc(); X <- X.update()
      W <- W.update()
      # print(sum(W))
      Hw <- lapply(1:N, H.update); tmp <- gc()
      
      obj <- obj.update();tmp <- gc()
      delta <- abs(obj0 - obj)/(mean(obj0, obj))
      obj0 <- obj
      objA[iters]<-obj
      plot(objA,main = paste0("Iteration no. ",iters))
      iters <- iters + 1
      setTxtProgressBar(pb, iters)
      print(delta)
      par(mfrow=c(1,2))
    }
    
    objAA[1:length(objA),i]<-objA
    setTxtProgressBar(pb, max.iters)
    if (iters == max.iters) {
      print("Warning: failed to converge within the allowed number of iterations. \n            Re-running with a higher max.iters is recommended.")
    }
    if (obj < best_obj) {
      W_m <- W
      Hw_m <- Hw
      H_m <- H
      A_m <- A
      best_obj <- obj
      best_seed <- rand.seed + i - 1
    }
    end_time <- difftime(Sys.time(), start_time, units = "auto")
    run_stats[i, 1] <- as.double(end_time)
    run_stats[i, 2] <- iters
    cat("\nConverged in ", run_stats[i, 1], " ", units(end_time), 
        ", ", iters, " iterations.\n", sep = "")
    if (print.obj) {
      cat("Objective:", obj, "\n")
    }
  }
  cat("Best results with seed ", best_seed, ".\n", sep = "")
  object@H <- H_m
  f1 <- function(H1){
    H1 <- lapply(1:length(object@scale.data), function(i){
      rownames(H1[[i]]) <- rownames(object@scale.data[[i]]);
      names(H1[[i]]) <- names(object@raw.data)
      return(H1[[i]])
    })
    return(H1)
  }
  f2<- function(m){colnames(m)<-object@var.genes;return(m)}
  object@H <- lapply(object@H,f1)
  Hw <- f1(Hw)
  W_m <- f2(W_m)
  object@W <- W_m
  colnames(object@W) = object@var.genes
  A <-lapply(A, f2)
  print("saving...")
  
  names(A)<-colnames(a)
  A<-lapply(names(A), function(x){
    A1<-A[[x]]
    rownames(A1)<-paste(x,1:nrow(A1),sep = "_")
    return(t(A1))})
  names(A)<-colnames(a)
  
  l<-list(A = A, W = W_m, H = H_m,Hw = Hw,
          obj = obj, k = k, k1 = k1,a = a,
          objAA = objAA, lambda.a = lambda.a,
          rand.seed = rand.seed, param = param)
  rownames(l$W)<-paste0("W",1:nrow(l$W))
  results<-list(A = A,W = t(l$W),param = l)
  
  if(missing(output.file)|is.null(output.file)){
    return(results)
  }
  saveRDS(results,file = output.file)
  return(results)
}

call_GMDF<-function(L,a,n.shared,n.spc,var.thresh = 0.3,output.file){
  # Input:
  # L (nx1) = list of n gene expression matrixes, one per dataset
  # a (nxk) = matrix of convariates, a(i,j) denoting covariate j in dataset i
  # n.shared (default 5): pre-defined number of shared programs
  # n.spc (default 5): pre-defined number of covariate-specific programs
  # var.thresh (default 0.3): Threshold used to identify variable genes.
  #                           Genes with expression variance greater than threshold (relative to mean) are selected.
  #                           (higher threshold -> fewer selected genes). Accepts single value or vector with separate var.thresh
  #                           for each dataset.
  # output.file (optional): The file where the results will be saved
  
  L1<-lapply(L,function(X) as(X,"sparseMatrix"))
  names(L1)<-names(L)
  
  ligerex = createLiger(L1) #Can also pass in more than 2 datasets
  ligerex = normalize(ligerex)
  ligerex = selectGenes(ligerex, var.thresh = var.thresh,combine = "intersection",do.plot = F)
  print(paste("Found",length(ligerex@var.genes),"var genes."))
  
  ligerex = scaleNotCenter(ligerex)
  print(length(ligerex@var.genes))
  print("Data size");print(dim(L[[1]]))
  print("Data size (liger)");print(dim(ligerex@raw.data[[1]]))
  print("Scaled data size");print(dim(t(ligerex@scale.data[[1]])))
  
  seed.id<-sample(1:1000,1)
  if(missing(output.file)){
    output.file <- NULL
  }else{
    output.file <- paste0(output.file,seed.id)
  }
  
  results<-GMDF(ligerex,k = n.shared,
                a = a,k1 = n.spc,
                print.obj = T,thresh = 1e-4,
                lambda.a = rep(1,ncol(a)),
                output.file = output.file,
                rand.seed = seed.id,nrep = 1)
  return(results)
}

solveNNLS <- function(C, B) {
  # liger::solve_nnls(C = C,B = B)
  # .Call('_liger_solveNNLS', PACKAGE = 'rliger', C, B)
  .Call('_rliger_solveNNLS', PACKAGE = 'rliger', C, B)
  
}

rbindlist<-function(l){
  X<-l[[1]]
  print(length(l))
  if(length(l)==1){
    return(X)
  }
  for(i in 2:length(l)){
    X<-rbind(X,l[[i]])
  }
  return(X)
}
