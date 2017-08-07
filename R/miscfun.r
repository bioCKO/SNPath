


Snp2Gene.Adist <- function(snp.info, gene.info, dist=20000){
    snp.info <- as.matrix(snp.info)
    gene.info <- as.matrix(gene.info)
    Chr.nam <- unique(as.vector(gene.info[,2]))
    gene.snplist <- list()
    gene.snplist[[nrow(gene.info)+1]] <- 0
    snps.pos <- as.numeric(snp.info[,3])
    
    start <- as.integer(as.vector(gene.info[,3]))
    end <- as.integer(as.vector(gene.info[,4])) 
    gene.pos <- cbind(start,end)
    gene.pos <- t(apply(gene.pos,1,range))
    gene.rang <- cbind(gene.pos[,1]-dist, gene.pos[,2]+dist)
       
    for (i in Chr.nam){
       idx.gene <- which(gene.info[,2]==i)
       idx.snp <- which(snp.info[,2]==i)
       snps.pos.i <- snps.pos[idx.snp]

       for (j in idx.gene){
           gene.rang.j <- gene.rang[j,]
           snps.j <- which(snps.pos.i>=gene.rang.j[1]&snps.pos.i<=gene.rang.j[2])
           if (length(snps.j)>0){
             gene.snplist[[j]] <- idx.snp[snps.j]
           }
       }
    }
    gene.snplist <- gene.snplist[1:nrow(gene.info)]
    names(gene.snplist) <- gene.info[,1]
    return(gene.snplist)  
}

gene2snp.Adist <- function(snp.info, gene.info, dist=20000){
    snp.info <- as.matrix(snp.info)
    gene.info <- as.matrix(gene.info)
    Chr.nam <- unique(as.vector(gene.info[,2]))
    snp.genelist <- list()
    snp.genelist[[nrow(snp.info)+1]] <- 0
    snps.pos <- as.numeric(snp.info[,3])
    
    start <- as.integer(as.vector(gene.info[,3]))
    end <- as.integer(as.vector(gene.info[,4])) 
    gene.pos <- cbind(start,end)
    gene.pos <- t(apply(gene.pos,1,range))
    gene.rang <- cbind(gene.pos[,1]-dist, gene.pos[,2]+dist)
       
    for (i in Chr.nam){
       idx.gene <- which(gene.info[,2]==i)
       idx.snp <- which(snp.info[,2]==i)
       snps.pos.i <- snps.pos[idx.snp]

       for (j in idx.gene){
           gene.rang.j <- gene.rang[j,]
           snps.j <- which(snps.pos.i>=gene.rang.j[1]&snps.pos.i<=gene.rang.j[2])
           if (length(snps.j)>0){
             oo <- idx.snp[snps.j]
             snp.genelist[oo] <- lapply(snp.genelist[oo], function(x,j){ x<-c(x, j)}, j=j)
           }
       }
    }
    snp.genelist <- snp.genelist[1:nrow(snp.info)]
    names(snp.genelist) <- snp.info[,1]
    return(snp.genelist)  
}


Snp2Gene.Rdist <- function(snp.info, gene.info, k=1){
    snp.info <- as.matrix(snp.info)
    gene.info <- as.matrix(gene.info)
    Chr.nam <- unique(as.vector(gene.info[,2]))
    gene.snplist <- list()
    gene.snplist[[nrow(gene.info)+1]] <- 0
    snps.pos <- as.numeric(snp.info[,3])
    for (i in Chr.nam){
       idx.gene <- which(gene.info[,2]==i)
       idx.snp <- which(snp.info[,2]==i)
       start <- as.integer(as.vector(gene.info[idx.gene,3]))
       end <- as.integer(as.vector(gene.info[idx.gene,4])) 
       mid <- (start+end)/2 
       od <- order(mid)
       for (j in idx.snp){
           snps.pos.i <- snps.pos[j]
           ooi <- which.min(abs(snps.pos.i-mid))
           odi <- which(od==ooi)+(-k:k)
           odi <- odi[odi>0]
           odi <- od[odi]
           gene_i_snps <- idx.gene[odi]
           gene_i_snps <- gene_i_snps[!is.na(gene_i_snps)]
           for (t in gene_i_snps){
             gene.snplist[[t]] <- c(gene.snplist[[t]], j)
           }
       }
    }
    gene.snplist <- gene.snplist[1:nrow(gene.info)]
    names(gene.snplist) <- gene.info[,1]
    return(gene.snplist)  
}


gene2snp.Rdist <- function(snp.info, gene.info, k=1){
    snp.info <- as.matrix(snp.info)
    gene.info <- as.matrix(gene.info)
    Chr.nam <- unique(as.vector(gene.info[,2]))
    snp.genelist <- list()
    snp.genelist[[nrow(snp.info)+1]] <- 0
    snps.pos <- as.numeric(snp.info[,3])
     
    for (i in Chr.nam){
       idx.gene <- which(gene.info[,2]==i)
       idx.snp <- which(snp.info[,2]==i)
       snps.pos.i <- snps.pos[idx.snp]
       
       start <- as.integer(as.vector(gene.info[idx.gene,3]))
       end <- as.integer(as.vector(gene.info[idx.gene,4])) 
       mid <- (start+end)/2
       od <- order(mid)
       
       for (j in idx.snp){
           snps.pos.i <- snps.pos[j]
           ooi <- which.min(abs(snps.pos.i-mid))
           odi <- which(od==ooi)+(-k:k)
           odi <- odi[odi>0]
           odi <- od[odi]
           gene_i_snps <- idx.gene[odi]
           gene_i_snps <- gene_i_snps[!is.na(gene_i_snps)]
           snp.genelist[[j]] <- gene_i_snps
       }
    }
    snp.genelist <- snp.genelist[1:nrow(snp.info)]
    names(snp.genelist) <- snp.info[,1]
    return(snp.genelist)  
}


Gene2Set <- function(gene, set){
    set.idx <- lapply(set, function(x){
        idx <- NULL    
        for (i in x) idx <- c(idx, which(gene==i))
        return(idx)
        })
    return(set.idx)
}



calc.fun <- function(cl=NULL, snp.dat, y, weights=NULL, snp.method=c("logiReg","chiSq")){
    
    ff=function(x,y, weights){ 
        if (is.null(weights)){
            mod <- glm(y~x, binomial(link = "logit"))
        } else{
            library(Zelig)
	    library(survey)
            mod <- zelig(cbind(y,1-y)~x,model="logit.survey", 
                           weights=~weights, data=data.frame(y,weights, x), cite=FALSE)    
        } 
        stat <- summary(mod)$coefficients[2,3:4]
        return(abs(stat))
    }

    ffc=function(x,y){ 
        mod <- chisq.test(x,y)
        return(c(mod$statistic, mod$p.value))
    }
    
    if (!is.null(cl)){
      if (snp.method[1]=="logiReg"|snp.method[1]=="l"){
        stat <- parRapply(cl, snp.dat, ff ,y=y, weights=weights)
      } else if (snp.method[1]=="chiSq"|snp.method[1]=="c"){
        stat <- parRapply(cl, snp.dat, ffc ,y=y)
      }
    } else {
      if (snp.method[1]=="logiReg"|snp.method[1]=="l"){
        stat <- apply(snp.dat, 1, ff ,y=y, weights=weights)
      } else if (snp.method[1]=="chiSq"|snp.method[1]=="c"){
        stat <- apply(snp.dat, 1, ffc ,y=y)
      }
    }
    stat <- matrix(stat, byrow=TRUE,ncol=2)
    return(list(pval=stat[,2], stat=stat[,1]))
}

