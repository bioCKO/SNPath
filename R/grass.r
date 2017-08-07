
#library(corpcor)

logistic.ridge.c <- function(X, y, grp, weights=NULL, lambda=1, tol=0.001, maxiter=100, intercept=TRUE, eq.eSNP=TRUE) { 
  grp <- as.integer(grp)
  grp <- grp-min(grp)+1
  if (intercept) {
      X <- cbind(1, X)
      grp <- c(1, grp+1)
  }
  n <- nrow(X)
  p <- ncol(X)
  ug <- unique(grp)
  ngrp <- length(ug)
  X.vector <- as.vector(t(X))
  Y.vector <- rep(1, length(y))
  uy <- sort(unique(as.numeric(y)))
  Y.vector[as.numeric(y)==uy[1]] <-  -1
  llik_output<- 0
  beta_output<-rep(0, p)
  intercept <- as.integer(intercept)
  eq.eSNP <- as.integer(eq.eSNP)
  
  if (is.null(weights) ) weights <- rep(1, length(y)) 
  
  results<-.C("logistic_ridge_c",
              as.double(X.vector),
              as.double(Y.vector),
              as.double(weights),
              as.integer(grp),
              as.integer(ngrp),
              as.integer(ug),
              as.integer(n),
              as.integer(p),
              as.double(lambda),
              as.double(tol),
              as.integer(maxiter),
              as.integer(intercept),
              as.integer(eq.eSNP),
              beta_output=as.double(beta_output),
              llik_output=as.double(llik_output),
              PACKAGE="SNPath"
              )
  list(beta=results$beta_output, ridge.llik=results$llik_output)
}

compute.llik.c <- function(X, y, weights, beta, intercept=TRUE) {
  if (intercept) {
      X <- cbind(1, X)
  }
  n <- nrow(X)
  p <- ncol(X)
  X.vector <- as.vector(t(X))
  Y.vector <- rep(1, length(y))
  uy <- sort(unique(as.numeric(y)))
  Y.vector[as.numeric(y)==uy[1]] <-  -1
  llik_output<- 0
 
  results<-.C("compute_loglik_c",
              as.double(X.vector),
              as.double(Y.vector),
              as.double(weights),
              as.double(beta),
              as.integer(n),
              as.integer(p),
              llik_output=as.double(llik_output),
              PACKAGE="SNPath"
              )
  return(results$llik_output)
}




find.lambda <- function(X, y, grp, weights=NULL, upper=1000, intercept=TRUE, eq.eSNP=TRUE, maxiter=100, tol= 0.001, method=c("AIC", "BIC")){
    method <- method[1]
    if (method=="AIC") {
        k <- 2 
    } else if (method=="BIC") {
        k <- log(length(y))
    }
    ABIC.fun <- function(lambda) {
      mod0 <- logistic.ridge.c(X, y, grp, weights=weights, lambda=lambda, intercept=intercept, eq.eSNP=eq.eSNP, tol=tol, maxiter=maxiter)
      aic <- -2*compute.llik.c(X, y, weights, mod0$beta, intercept=intercept)  + k*sum(mod0$beta!=0)
      return(aic)
    }

    lamr= 1:10
    lidx <- which.min(unlist(lapply(lamr, ABIC.fun),use.name=FALSE))
    
    while(lidx==10){
       lamr= lamr*10
       aic.val= unlist(lapply(lamr, ABIC.fun),use.name=FALSE)
       lidx <- which.min(aic.val)
    }
    lambda=lamr[lidx]
     
    return(lambda) 
}




gen.esnps <- function(snp.dat, snp.info, gene.info, gene.set, y, weights=NULL, gene.def=c("abs", "rel"),
                  k=1, dist=20000){
                  
    if (gene.def[1]=="abs"|gene.def[1]=="a"){
        gene.snplist <- Snp2Gene.Adist(snp.info, gene.info, dist=dist)
    } else if (gene.def[1]=="rel"|gene.def[1]=="r"){
        gene.snplist <- Snp2Gene.Rdist(snp.info, gene.info, k=k)
    }
    set.idx <- Gene2Set(gene.info[,1], gene.set)
    
    esnp.list <- list()
    gnam <- gene.info[,1]
    ug <- unique(unlist(set.idx))
    for (i in ug){
        z <- gene.snplist[[i]]
        if (is.null(z[1])){
          xvec <- NA
        } else if (length(z)>=3){
            ss <- fast.svd(snp.dat[z,])
            pi <- ss$d^2/sum(ss$d^2)  
            nsv <- max(sum(pi*length(z)>=1), 1)  
            xvec <- matrix(ss$v[, 1:nsv], ncol=nsv)
        } else {
            xvec <- t(snp.dat[z,])
            if (nrow(xvec)==1) xvec <- t(xvec)
        } 
        esnp.list[[i]] <- xvec
    }
    names(esnp.list) <- gene.info[1:length(esnp.list),1]
    
    
    if (is.null(weights)) weights <- rep(1,length(y))
    
    lis <- list()
    for (i in 1:length(gene.set)){
        oo <- unlist(esnp.list[set.idx[[i]]])
        oo <- oo[!is.na(oo)]
        dat <- matrix(oo,nrow=ncol(snp.dat))
        ng <- unlist(sapply(esnp.list[set.idx[[i]]], ncol))
        path.esnp <- dat
        gg <- rep(1:length(ng), ng)
        names(gg) <- rep(names(ng),ng)

        lis[[i]] <- list(path.esnp=dat, path.grp=gg, y=y, weights=weights)
    }
    names(lis) <- names(gene.set)

    return(lis)
}





getB.fun<-function(gseaLis, y0=NULL, weights0=NULL, lambda=NULL, nmsiz=100, intercept=TRUE, eq.eSNP=TRUE,tol=0.001,maxiter=100, method="AIC", return.lambda=TRUE){

     X=gseaLis$path.esnp
     grp=gseaLis$path.grp
     y=gseaLis$y
     weights= gseaLis$weights
     
     if (!is.null(y0)) {
       y<- y0
       weights <- weights0
     }
     
     ng=max(grp)
     nsplit= as.integer(ng/nmsiz)+1
     if (ng%%nmsiz==0) nsplit=nsplit-1
          
     sidx= split(1:ng, 1:nsplit)
     spidx=  lapply(sidx, function(x) {
         gidx=NULL
         for (j in x)
             gidx=c(gidx, which(grp==j))
             return(gidx)
     })
    
     if (is.null(lambda)) lamb=rep(NA, nsplit)
     beta=list()
     for (i in 1:nsplit){
         Xi=X[, spidx[[i]]] 
         grpi=grp[spidx[[i]]]
         ug=unique(grpi)
    
         if (is.null(lambda)){
             lamb[i] <- find.lambda(X=Xi, y=y,  grp=grpi, weights=weights,intercept=intercept,  
                     eq.eSNP=eq.eSNP,  method=method, tol=tol, maxiter=maxiter)
         } else {
             lamb=lambda
         }
             
         modi <- logistic.ridge.c(X=Xi, y=y, grp= grpi, weights=weights, lambda=lamb[i], intercept=intercept,
                    eq.eSNP=eq.eSNP, tol=tol,maxiter=maxiter)
         b=modi$beta[-1]
     
         for (u in ug){
                 beta[[u]]=b[which(grpi==u)]
         }
     }
     if (return.lambda){
       return(list(beta=beta, lambda=lamb))
     } else {
       return(beta)
     }
}

get0.fun<-function(seed, gseaLis, lambda, nmsiz=100, intercept=TRUE, eq.eSNP=TRUE,tol=0.001,maxiter=100, method="AIC", return.lambda=FALSE){
    set.seed(seed)
    y <- gseaLis$y
    weights <- gseaLis$weights
    idx=sample(1:length(y))
    y0 <- y[idx]
    weights0 <- weights[idx]
    
    beta <-  getB.fun(gseaLis, y=y0, weights0=weights0, lambda=lambda,intercept=intercept,  
                     eq.eSNP=eq.eSNP,  method=method, tol=tol, maxiter=maxiter, return.lambda=return.lambda)
    return(beta)
}


beta2stat<-function(beta){
  sapply(beta,function(x) sqrt(sum(x^2)))
}

getp<- function(stat,stat0, nominal=FALSE){
   m=rowMeans(stat0)
   var=apply(stat0,1,var)
   s=sqrt(sum((stat-m)^2/var))
   s0=apply(stat0,1,function(x) (x-mean(x))^2/var(x))
   s0=sqrt(rowSums(s0))
   
   if (nominal==FALSE) {
     p=mean(s0>=s)
   } else {
     p=1-pnorm(s,mean(s0),sd(s0))
   }
   return(p)
}

## require load the library corpcor and SNPath to both local and all clusters

grass <- function(cl=NULL, snp.dat, snp.info, gene.info, gene.set, y, weights=NULL, 
                 gene.def=c("abs", "rel"), dist=20000,  k=1, B=100, nominal.p=FALSE){
 

  gsea.obj <- gen.esnps(snp.dat=snp.dat, snp.info=snp.info, gene.info=gene.info, 
                   gene.set=gene.set, y=y, weights=weights, gene.def=gene.def[1], k=k, dist=dist)
  
  if (is.null(cl)){
      res<-lapply(gsea.obj, getB.fun)
  } else {
      res<-parLapply(cl, gsea.obj, getB.fun)
  }
  lambda=lapply(res, function(x) x$lambda)
  betas= lapply(res, function(x) x$beta)

  beta0s <- list()
  Np <- length(gsea.obj)
  for (i in 1:Np){
      if (is.null(cl)){
          beta0s[[i]] <- lapply(1:B, get0.fun, gseaLis=gsea.obj[[i]], lambda=lambda[[i]], return.lambda=FALSE)
      } else {  
          beta0s[[i]] <- parLapply(cl, 1:B, get0.fun, gseaLis=gsea.obj[[i]], lambda=lambda[[i]], return.lambda=FALSE)
      }
  }

  stats = lapply(betas, beta2stat)
  stat0s = lapply(beta0s, function(x) matrix(sapply(x, beta2stat), byrow=F, ncol=B))

  pval <- rep(NA, length(stats))
  names(pval) <- names(stats)
  for (i in 1:length(stats)){
    pval[i] <- getp(stats[[i]], stat0s[[i]],nominal=nominal.p)
  }

  return(pval)
}







