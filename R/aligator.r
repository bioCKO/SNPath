




aligator <-  function(cl=NULL, snp.info, gene.info, gene.set, snp.pval,
      gene.def=c("abs", "rel"), dist=20000, k=1, Bresample=5000, snp.pcut=0.005){


    if (gene.def[1]=="abs"|gene.def[1]=="a"){
        snp.genelist <- gene2snp.Adist(snp.info, gene.info, dist=dist)
    } else if (gene.def[1]=="rel"|gene.def[1]=="r"){
        snp.genelist <- gene2snp.Rdist(snp.info, gene.info, k=k)
    }
    snp.genelist <- lapply(snp.genelist, as.vector)
    names(snp.genelist) <- snp.info[,1]    
    set.idx <- Gene2Set(gene.info[,1], gene.set)
    
    idx <- which(snp.pval<=snp.pcut)
    hit.gene <- unique(unlist(snp.genelist[idx], use.name=FALSE))   ## this is the top hit genes from observed data at the cutoff
    pi <- length(idx)
    N <- length(hit.gene)
 
    p <- length(snp.pval)

    fun.aligator <- function(seed, snp.genelist, p, pi, N){
        set.seed(seed)
        idxi <- sample(1:p)
        hit.gene0 <- unique(unlist(snp.genelist[idxi[1:pi]], use.name=FALSE))
        N0 <- length(hit.gene0)
        while(N0<N){
            pi <- pi+1            
            hit.gene0 <- unique(unlist(snp.genelist[idxi[1:pi]], use.name=FALSE))
            N0 <- length(hit.gene0)
        }
        while(N0>N){
	    pi <- pi-1            
	    hit.gene0 <- unique(unlist(snp.genelist[idxi[1:pi]], use.name=FALSE))
	    N0 <- length(hit.gene0)
        }
        return(hit.gene0)
    }


    if (is.null(cl)){
        hit.gene0 <- lapply(1:Bresample, fun.aligator, snp.genelist=snp.genelist, p=p, pi=pi, N=N)  ## those are the hitted gene indices from the re-sampled data
    } else {
        hit.gene0 <- parLapply(cl, 1:Bresample, fun.aligator, snp.genelist=snp.genelist, p=p, pi=pi, N=N)
    }
    
    pval <- NULL
    for (i in 1:length(set.idx)){
        path <- set.idx[[i]]
        count0 <- sapply(hit.gene0, function(hit, path){
            sum(table(c(path,hit))>1)    ## I "table" the indice of the genes in the pathway and the indice of genes from top hits and see how many overlaps
        }, path=path)
        count <- sum(table(c(path, hit.gene))>1)
        pval <- c(pval, mean(count0>=count))
    }
    names(pval) <- names(set.idx)

    return(pval)
}


