


\name{ aligator }
\alias{ aligator }
\title{ aligator }

\description{
    This function calculates the p-value for the association between a pathway
    and disease risk using the ALIGATOR algorithm.
}

\usage{
    aligator(cl=NULL, snp.info, gene.info, gene.set, snp.pval,
      gene.def=c("abs", "rel"), dist=20000, k=1, Bresample=5000, snp.pcut=0.005)
}

\arguments{
  \item{ cl }{Cluster object for parallel computing. We recommand using the
    `snow' and `Rmpi' libraries to generate the clusters. See examples
    for illustration. If clusters are not provided, cl=NULL, computing will
    be conducted on the local machine. Parallel computing is optional but highly recommended. 
  }
  
  \item{ snp.info }{A matrix with three columns: snp.Name (name of SNP), chr (chromosome),
    and pos (position where the SNP is located). It is suggested that the SNPs are
    sorted by their locations in the genome. Please insure that the same reference is used
    for SNP position and for gene position.
  }

  \item{ gene.info }{A matrix with four columns with the order of gene.Name, chr,
   Start, and End. The name of the gene is under gene.Name. The
   chromosome is chr, and the start and end coordinate (start < end) are
   in the Start and End columns. It is suggested that the genes are
    sorted by their locations in the genome.
  }
  
    
  \item{ gene.set }{A list of pathways. Each element of the list consists of the gene
    names of a pathway of interest.
  }

  \item{ snp.pval }{A numerical vector of individual SNP association p-values. Aligator
    takes the SNP p-values as input and does not need individual level SNP data. If
    the snp data is available, one can also use \code{calc.fun} to calculate p-values and
    associated statistics.
  }

  \item{ gene.def }{An optional character string giving a method for assigning SNPs to 
    the nearest genes by genomic distance. It can be "abs", which stands for absolute 
    genomic distance in terms of bp, or "rel" for relative distance. 
  }
  
  \item{ dist }{A numerical value of absolute genomic distance (in base pair). Only used when 
    \link{gene.def} = "abs". For example, gene.def="abs" and dist=20000, SNPs that are within
    20kb of either end of a gene will be assigned to that gene.
  }
  
  \item{ k }{A numerical value of relative genomic distance. Only used when \link{gene.def} =
    "rel". Each SNP is assigned to the k nearest genes on either side plus the
    gene that contains the SNP.  A SNP can be assigned to a maximum of \eqn{2k+1}{2*k+1} genes.
  }
  
  \item{ Bresample }{The number of resamplings performed to calculate null statistics. 
  }
  
  \item{ snp.pcut }{A numerical value. Only the SNPs with marginal p-value \code{snp.pval} less than
    this cutoff are used in the analysis.
  }


}

\seealso{
    \code{  \link{ gseaSnp } } , \code{  \link{ grass } }, \code{  \link{ plinkSet } } 
}

\details{
    This algorithm takes the individual SNP association p-values as input and uses a preselected
    p-value threshold to define a set of significantly associated SNPs. It then counts the number
    of genes in the pathway that contains these SNPs, with each gene counted only once, regardless
    of the number of significant SNPs in the gene. It uses resampling of SNPs to establish the
    null distribution. Thus it only requires a p-value or summary statistic from each SNP as input,
    and can be used when individual level SNP data are not available.

}


\value{
  Returns a numerical vector of p-values, corresponding to the pathways in \link{gene.set}.  
}

\references{
  Holmans et al. (2010). Gene ontology analysis of GWA study data sets provides insights into
  the biology of bipolar disorder. American Journal of Human Genetics, 85(1): 13-24.
}


\examples{
  data(simDat)
  
  library(snow)
  cl <- makeCluster(c("localhost","localhost"), type = "SOCK")

  pval <- calc.fun(cl=cl, snp.dat=snpDat, y=y, snp.method="logiReg")$pval

  path.pval <- aligator(cl=cl, snp.info=snp.info, gene.info=gene.info,
                   gene.set=sim.pathway, snp.pval=pval, gene.def="rel", snp.pcut=0.05)
 
  stopCluster(cl)
  
}

%\keyword{ TBD }

