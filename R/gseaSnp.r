


######################################################################################################################
######################################################################################################################
## the following two functions are the ones for running Wang et al. (2007) method.

## this function calculate the Enrichment statistics for a each of the pathways
## path.idx is a list, each of its element is the indice of a set of genes in that pathway
## gene.stats is a vector of the chi-square statistics from all SNPs.
## gene.snplist is a list for all genes, each element is the indice of snps for that gene. It is the result returned by
## the function Snp2gene

  getES <- function(path.idx, gene.stats, gene.snplist, p=1){
    ES=rep(0,length(path.idx)) 
    rk <- rank(-gene.stats,ties.method="first")
    N=length(gene.snplist)  ## total number of genes

    for (i in 1:length(path.idx)){
      set.idx=path.idx[[i]]  ## the gene indice for this i-th pathway
      Nh=length(set.idx)  ## number of genes in this pathway
      oo <- sort(rk[set.idx])
  
      ES.all= -(1:N)/(N-Nh)
      statj=gene.stats[set.idx]
      statj=-sort(-statj)
      Nr=sum(abs(statj)^p)
      for (j in 1:(Nh-1)){
         jj=sum(abs(statj[1:j])^p)
         ES.all[oo[j]:(oo[j+1]-1)] = ES.all[oo[j]:(oo[j+1]-1)]+jj/Nr+j/(N-Nh)
      }   
      ES.all[N]=0
      ES[i]=max(ES.all)
    }
    return(ES)
  }

getp.gsea <- function(path.idx, stats, stats0, gene.snplist){
  gene.stats=sapply(gene.snplist, function(x, stat){
    max(stat[x])
   },stat=stats)
  ES=getES(path.idx=path.idx,gene.stats=gene.stats, gene.snplist=gene.snplist)

  ES0=NULL
  for (i in 1:ncol(stats0)){
    s0=sapply(gene.snplist, function(x, stat){
       max(stat[x])
       },stat=stats0[,i])
     ES0=cbind(ES0,getES(path.idx,gene.stats=s0, gene.snplist=gene.snplist))
  }

  mES=apply(ES0,1,mean)
  sdES=apply(ES0,1,sd)
  NES=(ES-mES)/sdES

  NES0=apply(ES0,2, function(x,mES,sdES) (x-mES)/sdES, mES=mES,sdES=sdES)

  pval=rep(0,length(NES))
  for (i in 1:length(NES)){
    pval[i]=mean(NES0>=NES[i])
  }

  #FDR=rep(0,length(NES))
  #for (i in 1:length(NES)){
  #  FDR[i]=mean(NES0>=NES[i])/mean(NES>=NES[i])
  #}
  #FDR[FDR>1]=1
  #names(pval) <- names(FDR) <- names(path.idx)
  
  return(pval)
}






gseaSnp <- function(cl=NULL, snp.dat, snp.info, gene.info, gene.set, y, weights=NULL, 
       snp.method=c("logiReg","chiSq"), gene.def=c("abs", "rel"), dist=20000, k=1, B=100){


    if (gene.def[1]=="abs"|gene.def[1]=="a"){
        gene.snplist <- Snp2Gene.Adist(snp.info, gene.info, dist=dist)
    } else if (gene.def[1]=="rel"|gene.def[1]=="r"){
        gene.snplist <- Snp2Gene.Rdist(snp.info, gene.info, k=k)
    }
    gene.snplist <- lapply(gene.snplist, as.vector)
    names(gene.snplist)=gene.info[,1]
    gene.snplist <- gene.snplist[sapply(gene.snplist,length)>=1]  ## we only keep the genes have at least 1 SNP

    set.idx <- Gene2Set(names(gene.snplist), gene.set)

    ii <- unique(unlist(set.idx,use.name=F))   ## we only calculate statistics for the snps that being used in any of the pathways
    snp.idxall=unique(unlist(gene.snplist[ii],use.name=F))

    if (!is.null(weights)){
       clusterEvalQ(cl, {
	     library(Zelig)
	     library(survey)
       })
    }
    
    ff=function(x,y, weights){ 
        if (is.null(weights)){
            mod <- glm(y~x, binomial(link = "logit"))
        } else{
            mod <- zelig(cbind(y,1-y)~x,model="logit.survey", weights=~weights, data=data.frame(y,weights, x), cite=FALSE)    
        } 
        stat <- summary(mod)$coefficients[2,3]
        return(abs(stat))
    }

    ffc=function(x,y){ 
            stat <- chisq.test(x,y)$statistic
            return(abs(stat))
    }
    
    stats=rep(0,nrow(snp.dat))
    dat <- snp.dat[snp.idxall, ]
    if (is.null(cl)){
        if (snp.method[1]=="logiReg"|snp.method[1]=="l"){
            stats[snp.idxall] <- apply(dat, 1, ff ,y=y, weights=weights)
        } else if (snp.method[1]=="chiSq"|snp.method[1]=="c"){
            stats[snp.idxall] <- apply(dat, 1, ffc ,y=y)
        }
    } else {
        if (snp.method[1]=="logiReg"|snp.method[1]=="l"){
            stats[snp.idxall] <- parRapply(cl, dat, ff ,y=y, weights=weights)
        } else if (snp.method[1]=="chiSq"|snp.method[1]=="c"){
            stats[snp.idxall] <- parRapply(cl, dat, ffc ,y=y)
        }
    }
    
    stats0=matrix(0,nrow=nrow(snp.dat),ncol=B)
    for (i in 1:B){
        set.seed(i)
        ii <- sample(1:length(y))
        y0 <- y[ii]
        weights0 <- weights[ii]

        if (is.null(cl)){
            if (snp.method[1]=="logiReg"){
        	stats0[snp.idxall, i] <- apply(dat, 1, ff ,y=y0, weights=weights0)
	    } else if (snp.method[1]=="chiSq"){
	        stats0[snp.idxall, i] <- apply(dat, 1, ffc ,y=y0)
            }
        } else {
            if (snp.method[1]=="logiReg"){
        	stats0[snp.idxall, i] <- parRapply(cl, dat, ff ,y=y0, weights=weights0)
	    } else if (snp.method[1]=="chiSq"){
	        stats0[snp.idxall, i] <- parRapply(cl, dat, ffc ,y=y0)
            }
        }
        
    }

    pval <- getp.gsea(set.idx, stats, stats0, gene.snplist) 
    names(pval) <- names(set.idx)
    return(pval)
}

