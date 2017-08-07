

statFun <- function(snpStat, snp.dat, snp.pcut=0.005, snp.r2cut=0.5){

  snp.pval <- snpStat$pval
  snp.stat <- snpStat$stat
  n <- ncol(snp.dat)
  ii <- which(snp.pval<=snp.pcut) 
  if (length(ii)>0){
    dat3 <- matrix(snp.dat[ii,], byrow=T, ncol=n)
    stat3 <- snp.stat[ii]
    statR <- NULL
    while(length(dat3)>0){
      od <- order(-stat3)
      x1 <- matrix(dat3[od[1],], byrow=T, ncol=n)
      ni <- length(od)
    
      statR <- c(statR, stat3[od[1]])
    
      idx2 <- od[2:ni]
      x2 <- matrix(dat3[idx2, ], byrow=T, ncol=n)
    
      corr <- apply(t(x2), 2, cor, y=t(x1), use="pairwise")
      oo <-  which(abs(corr)>=snp.r2cut)
      if (length(oo)>0) {
        dat3 <- matrix(dat3[-idx2[oo],], byrow=T, ncol=n)
        stat3 <- stat3[-idx2[oo]]
      }

      dat3 <- matrix(dat3[-od[1],],ncol=n)
      stat3 <- stat3[-od[1]]
    }  
    fstat <- mean(statR)
  } else {
    fstat <- 0
  }
  return(fstat)
  
}




plinkSet <- function(cl=NULL, snp.dat, snp.info, gene.info, gene.set, y, weights=NULL, 
                       snp.method=c("logiReg","chiSq"),gene.def=c("abs", "rel"), 
                       dist=20000,  k=1, B=100, snp.pcut=0.005, snp.r2cut=0.5){

    if (gene.def[1]=="abs"|gene.def[1]=="a"){
        gene.snplist <- Snp2Gene.Adist(snp.info, gene.info, dist=dist)
    } else if (gene.def[1]=="rel"|gene.def[1]=="r"){
        gene.snplist <- Snp2Gene.Rdist(snp.info, gene.info, k=k)
    }
    gene.snplist <- lapply(gene.snplist, as.vector)
    names(gene.snplist)=gene.info[,1]
    gene.snplist <- gene.snplist[sapply(gene.snplist,length)>=1]  ## we only keep the genes have at least 1 SNP

    set.idx <- Gene2Set(names(gene.snplist), gene.set)

    n <- length(y)
    nsnp <- nrow(snp.dat)

    vv <- unique(unlist(set.idx, use.name=F))
    snp.idx <- unlist(gene.snplist[vv], use.name=F)
    snp.idx <- unique(snp.idx)
        
    oo <- calc.fun(cl=cl, snp.dat=snp.dat[snp.idx,], y=y, weights=weights, snp.method=snp.method)
    pval <- rep(0, nsnp); pval[snp.idx] <- oo$pval
    stat <- rep(0, nsnp); stat[snp.idx] <- oo$stat

    path.snpidx <- lapply(set.idx, function(ii, gene.snplist){
      sidx <- unlist(gene.snplist[ii], use.name=F)
      sidx <- unique(sidx)
      return(sidx)
    }, gene.snplist=gene.snplist)
    

    permlist <- list()
    permlist[[B]] <- 0
    for (b in 1:B){
       set.seed(b)
       perm <- sample(1:n)
       y0 <- y[perm]
       weights0 <- weights[perm]
       oo <- calc.fun(cl=cl, snp.dat=snp.dat[snp.idx,], y=y0, weights=weights0, snp.method=snp.method)
       pval0 <- rep(0, nsnp); pval0[snp.idx] <- oo$pval
       stat0 <- rep(0, nsnp); stat0[snp.idx] <- oo$stat
       permlist[[b]] <- list(pval=pval0,stat=stat0)
    }

    path.pval <- NULL
    for (i in 1:length(set.idx)){
        snpi <- path.snpidx[[i]]
        dat <- snp.dat[snpi, ]
        snpStat <- list(pval=pval[snpi], stat=stat[snpi])
        
        pstat <- statFun(snpStat, snp.dat=dat, 
                     snp.pcut=snp.pcut, snp.r2cut=snp.r2cut)
    
        snpStat.lis <- lapply(permlist, function(x, snpi){
           return(list(pval=x$pval[snpi],stat=x$stat[snpi]))   
        }, snpi=snpi)
        
        if (is.null(cl)){
            pstat0 <- lapply(snpStat.lis, statFun,  snp.dat=dat,
                    snp.pcut=snp.pcut, snp.r2cut=snp.r2cut)  
        } else {
            pstat0 <- parLapply(cl, snpStat.lis, statFun,  snp.dat=dat,
                    snp.pcut=snp.pcut, snp.r2cut=snp.r2cut) 
        }
        pstat0 <- unlist(pstat0, use.name=F)
    
        path.pval <- c(path.pval, mean(pstat0>=pstat))
    }
    names(path.pval) <- names(gene.set)
    
    return(path.pval)
}