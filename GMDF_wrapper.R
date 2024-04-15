GMDF_wrapper<-function(E,a, k = 5, k1 = 3,N1 = 1,outputdir,
                       lambda.a = rep(1,ncol(a)),var.thresh = 0.2,combine = "intersection"){
  # Annotations:
  #             g = number of genes
  #             k = number of shared programs
  #             k1 = number of context-specific programs
  #             n = number of datasets
  #             m = number of contexts (e.g, cancer type)
  
  # Input:
  # E: list of n gene expression matrices, one per dataset/condition. 
  # a: n x m matrix denoting the value of m contexts in n datasets provided in E.
  # k: Initial estimate of the number of shared programs.
  # k1: Initial estimate of number of context-specific programs per context.
  # N1: Number of times to run GMDF. If N1 > 1 then multiple solutions will be combined.
  # outputdir: The output directory to save the single GMDF run results. Required if N1 > 1.
  
  # Output when N1 > 1 are the results of N1 combined GMDF runs:
  #     - sig – the final shared signatures (top 100 genes in the programs represented in the Wf matrix).
  #     - Wf – Wf – the final GMDF shared programs based on all the GMDF runs, such that Wf[i,j] denotes the weight of gene i in program j.
  #     - W.multi.run – weight of genes (rows) in each of the shared programs identified in the single GMDF runs (columns)
  #     - W.clusters – the clustering annotation of each one of the shared GMDF programs from the single runs
  
  # Output description for N1 = 1, the results of a single GMDF run:
  # The decomposition of the data
  # Et[[i]] ~ (Hw[[i]]) x W + ∑j(a[i,j] x  H[[j]][[i]] x A[[j]])
  # Where Et is transpose(E[[i]][g,])
  # For additional information about the output and its interpretation
  # see Figure 1A and equation (1) as provided in the STAR METHODS describing GMDF.
  
  # Output when N1 = 1:
  # Hw - shared programs usage per dataset: list of n matrix, one per dataset, each of size g x k
  # W - g x k matrix representing the shared programs.
  # A - context-specific programs, k1 x g matrix per each of the m contexts
  # H[[j]] - context-specific program usage per context j, with a (g x k1) matrix per each of the n datasets
  # obj - the final value of the objective function, minimizing the reconstruction error
  # objAA - the values of the objective function in each iteration; NA in case the run terminated before the 100th iteration.
  # param - the parameters that were provided as input and used
  # input.data - the input object
  
  seed.id<-sample(1:1000000,1)
  set.seed(seed.id)
  
  ligerex = createLiger(E)
  ligerex = normalize(ligerex)
  ligerex = selectGenes(ligerex, var.thresh = var.thresh,combine = combine,do.plot = F)
  ligerex = scaleNotCenter(ligerex)
  
  if(N1==1){
    rslts<-GMDF(ligerex,k = k,a = a,k1 = k1,
                print.obj = T,thresh = 1e-4,lambda.a = lambda.a,
                rand.seed = seed.id,nrep = 1)
    return(rslts)
  }
  
  rslts<-lapply(1:N1, function(i){
    seed.id<-sample(1:1000000,1)
    set.seed(seed.id)
    rslts<-GMDF(ligerex,k = k,a = a,k1 = k1,
                print.obj = T,thresh = 1e-4,lambda.a = lambda.a,
                rand.seed = seed.id,nrep = 1,saveFile = paste0(outputdir,"/GMDF_run",i,".rds"))
    return(rslts)
  })
  
  W<-lapply(rslts, function(x) return(x$W))
  rslts1<-GMDF_combine(W)
  return(rslts1)
}

GMDF<-function (object, k, k1, a, lambda.a, thresh = 1e-04, max.iters = 100, 
                nrep = 1, A.init = NULL,H.init = NULL, W.init = NULL, V.init = NULL,
                rand.seed = 1, print.obj = F,int.val = c(0.1,0.1),saveFile){
  
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
      Hw <- lapply(1:N, H.update); tmp <- gc()
      
      obj <- obj.update();tmp <- gc()
      delta <- abs(obj0 - obj)/(mean(obj0, obj))
      obj0 <- obj
      objA[iters]<-obj
      plot(objA,main = paste0("Iteration no. ",iters),xlab = "Iteration",ylab = "Objective value")
      iters <- iters + 1
      setTxtProgressBar(pb, iters)
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
      return(H1[[i]])
    })
    names(H1) <- names(object@raw.data)
    return(H1)
  }
  f2 <- function(m){
    colnames(m)<-object@var.genes
    rownames(m)<-paste0("A",1:nrow(m))
    return(m)
  }
  object@H <- lapply(object@H,f1)
  Hw <- f1(Hw)
  H_m<-lapply(H_m, f1)
  W_m <- f2(W_m)
  rownames(W_m)<-paste0("W",1:nrow(W_m))
  object@W <- W_m
  A <-lapply(A, f2)
  
  names(Hw)<-rownames(a) # Shared programs usage per dataset
  names(A)<-colnames(a) # Context-specific program usage per context, as a list with matrix per dataset
  names(H_m)<-colnames(a) # Context-specific program usage per context, as a list with matrix per dataset
  
  output<-list(A = A, W = t(W_m), H = H_m,Hw = Hw,
               obj = obj, objAA = objAA,
               param = list(k = k, k1 = k1,a = a,lambda.a = lambda.a,seed = rand.seed),
               input.data = object)
  
  # For interpretation of the output see Figure 1A and equation (1)
  # as provided in the STAR METHODS describing GMDF.
  if(!missing(saveFile)){
    saveRDS(output,file = saveFile)
  }
  return(output)
}

solveNNLS <- function(C, B) {
  # liger::solve_nnls(C = C,B = B)
  # .Call('_liger_solveNNLS', PACKAGE = 'rliger', C, B)
  .Call('_rliger_solveNNLS', PACKAGE = 'rliger', C, B)
  
}

rbindlist<-function(l){
  X<-l[[1]]
  if(length(l)==1){
    return(X)
  }
  for(i in 2:length(l)){
    X<-rbind(X,l[[i]])
  }
  return(X)
}

GMDF_combine_pancancer<-function(){
  print("Combining pan-cancer CD8 T cell GMDF shared programs.")
  rslts<-readRDS("GMDF_W_pancanCD8.rds")
  rslts1<-GMDF_combine(rslts)
  return(rslts1)
}

GMDF_combine<-function(rslts,path2output){
  
  if(missing(rslts)){
    files<-list.files(path2output,full.names = T)
    rslts<-lapply(files,function(x){x<-readRDS(x);return(x$W)})
  }
  
  names(rslts)<-paste0("R",1:length(rslts))
  genes<-table(unlist(lapply(rslts, function(W) rownames(W))))
  g<-names(genes)[genes>=(length(rslts)*0.95)]
  W<-NULL
  for(x in names(rslts)){
    W1<-rslts[[x]]
    colnames(W1)<-paste(x,1:ncol(W1),sep = "_")
    W1<-W1[match(g,rownames(W1)),]
    W<-cbind(W,W1)
  }
  W<-W[,colSums(is.na(W))<100]
  W<-W[rowSums(is.na(W))==0,]
  W.clusters<-mNMF_cluster(W,method = "euclidean")
  b<-get.abundant(W.clusters,3,boolean.flag = T)
  Wf<-apply(W,2,function(x) x/max(x))
  Wf<-t(average.mat.rows(t(Wf[,b]),W.clusters[b]))
  sig<-get.top.elements(-Wf)
  sig.w<-split(colnames(W[,b]),W.clusters[b])
  names(sig.w)<-colnames(Wf)
  
  sig.w<-sig.w[order(laply(sig.w,length),decreasing = T)]
  sig<-sig[names(sig.w)]
  
  over.sig<-GO.enrichment.lapply(sig,names(genes),sig,valuType = 3)
  diag(over.sig)<-0
  x<-1
  while(x<ncol(over.sig)){
    b<-over.sig[,x]<20
    over.sig<-over.sig[b,b]
    x<-x+1
  }
  idx<-colnames(over.sig)
  sig.w<-sig.w[idx]
  sig<-sig[idx]
  rslts1<-list(sig = sig,
               Wf = Wf,
               W.multi.run = W,
               W.clusters = W.clusters)
  return(rslts1)
  
}

mNMF_cluster<-function(W,res = 1.5,method = "euclidean",gene.ids,cexRow = 0.5,ylab = "Genes"){
  
  if(is.element(method,c("pearson","spearman"))){
    d.row<-as.dist(1-cor(t(W), method = method))
    d.col<-as.dist(1-cor(W, method = method))
  }else{
    d.row<-dist(W,method = method)
    d.col<-dist(t(W),method = method)
  }
  
  hr.row <- hclust(d.row, method="complete")
  hr.col <- hclust(d.col, method="complete")
  mycl.col <- cutree(hr.col, h=max(hr.col$height/res))
  mycl.col<-paste0("Wf",mycl.col)
  
  return(mycl.col)
}

get.abundant<-function(v,abn.c = 2,boolean.flag = F,top,decreasing = T){
  m<-as.matrix(table(v))
  m<-as.matrix(m[order(m,decreasing = decreasing),])
  if(!missing(top)){
    abn.c<-m[top]
  }
  m<-m[m>=abn.c,]
  abn.names<-names(m)
  if(boolean.flag){
    b<-is.element(v,abn.names)
    return(b)
  }
  return(abn.names)
}

average.mat.rows<-function(m,ids,f = colMeans){
  ids.u<-sort(unique(ids))
  m1<-get.mat(ids.u,colnames(m))
  for(x in ids.u){
    b<-is.element(ids,x)
    if(sum(b)==1){
      m1[x,]<-m[b,]
    }else{
      m1[x,]<-f(m[b,])
    }
  }
  return(m1)
}

get.mat<-function(m.rows,m.cols,data = NA){
  m<-matrix(data = data, nrow = length(m.rows),ncol = length(m.cols),
            dimnames = list(m.rows,m.cols))
  return(m)
}

get.top.elements<-function (m,q = 100,min.ci = NULL,main = "",sort.flag = T){
  top.l<-list()
  v<-rownames(m)
  for (i in 1:ncol(m)){
    mi<-m[,i];mi<-mi[!is.na(mi)]
    idx<-order(mi,decreasing = F)
    ci <- mi[idx[min(q,length(mi))]]
    ci <- min(ci,min.ci)
    b <- m[,i]<=ci
    b[is.na(m[,i])]<-F
    if(sort.flag){
      top.l[[i]]<-sort(v[b])
    }else{
      top.l[[i]]<-v[b][order(m[b,i])]
    }
    
  }
  if(main!=""){main<-paste0(main,".")}
  names(top.l)<-paste0(main,colnames(m))
  return(top.l)
}

GO.enrichment.lapply<-function(go.env,genes,sig,valuType = 1){
  m<-t(laply(sig,function(x) GO.enrichment(go.env,genes,x)[,valuType]))
  colnames(m)<-names(sig)
  return(m)
}

GO.enrichment<-function(go.env,genes,selected.genes,add.genes = F,sort.flag = F){
  b<-is.element(genes,selected.genes)
  p<-laply(go.env,function(x) get.hyper.p.value(is.element(genes,x),b))
  rownames(p)<-names(go.env)
  if(add.genes){
    p<-cbind.data.frame(p,
                        Genes = laply(go.env,function(x) paste(intersect(selected.genes,x),collapse = ", ")))
  }
  if(sort.flag){
    p<-p[order(p[,1]),]
  }
  return(p)
}

get.hyper.p.value<-function(b1,b2,full.flag = T){
  p1<-NA;p2<-NA;e<-0;
  if(any(b1)&&any(b2)){
    p1<-max(1-phyper(sum(b1&b2)-1, sum(b1), sum(!b1), sum(b2)),1e-17)
    e<-sum(b2)*(sum(b1)/length(b2))
    p2<-sum(b1&b2)/e
  }
  if (full.flag){
    p<-c(p1,p2,sum(b1&b2),e)
    names(p)<-c('hyper.p.value','ob.vs.exp','ob','exp')
    return(p)  
  }
  return(p1)
}
