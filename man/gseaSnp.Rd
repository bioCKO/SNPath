

\name{ gseaSnp }
\alias{ gseaSnp }
\title{ gseaSnp }

\description{
    This function calculates the p-value of disease-association for a pathway by
    the algorithm proposed in Wang et al. (2007) AJHG.
}

\usage{
    gseaSnp(cl=NULL, snp.dat, snp.info, gene.info, gene.set, y, weights=NULL, 
       snp.method=c("logiReg","chiSq"), gene.def=c("abs", "rel"), dist=20000, k=1, B=100)
}

\arguments{
  \item{ cl }{Cluster object for parallel computing. We recommand using the
    `snow' and `Rmpi' libraries to generate the clusters. See examples
    for illustration. If clusters are not provided, cl=NULL, computing will
    be conducted on the local machine. Parallel computing is optional but highly recommended. 
  }
  
  \item{ snp.dat }{A matrix with snp data. Data are coded as 0, 1, 2, corresponding
    to homozygotes for the major allele, heterozygotes, and homozygotes for 
    the minor allele, respectively. Each row of the matrix is one SNP and each column
    represents one subject.
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

  \item{ y }{A vector of case/control status. Cases are coded as 1 and controls are 
    coded as 0.
  }
  
  \item{ weights }{An optional numerical vector specifying each subject’s sample weight, the
    inverse probability that the subject is selected. Use of weights requires loading the 
    libraries `Zelig' and `survey'.
  }

  \item{ snp.method }{An optional character string which specifies the method used to compute
    statistics for the marginal association analysis of each individual SNP. Two options are provided.
    By default, snp.method ="logiReg", logistic regression (additive effect model) will 
    be used. If snp.method = "chiSq", chi-square test will be conducted.  
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
  
  \item{ B }{The number of permutations performed to calculate null statistics. 
  }
 

}

\seealso{
  \code{  \link{ aligator } } , \code{  \link{ grass } }, \code{  \link{ plinkSet } } 
}

\details{
  This algorithm works as follows:
 1. Assign SNPs to genes based on absolute or relative genome location.
 2. Use the top individual SNP association statistic within the gene as the statistic of the gene.
 3. Rank all the genes by significance.
 4. Compares the distribution of the ranks of genes from a given pathway to that of the remaining
    genes via a weighted Kolmogorov-Smirnov test. Greater weights are given to genes with more 
    extreme statistic values.
 5. Permute the case/control status B times and repeat the above procedure to obtain null pathway
    statistics.
 6. The p-value for the pathway is calculated as the frequency of null statistics exceeding the 
    observed ones.
 
}

\value{
  Returns a numerical vector of p-values, corresponding to the pathways in \link{gene.set}.  
}

\references{
  Kai Wang, Mingyao Li and Maja Bucan (2007). Pathway-based approaches for analysis of 
  genomewide association studies. American Journal of Human Genetics, 81(6): 1278-1283.
}

%\author{
%  TBD
%}

%\note{
%  TBD
%}

\examples{
  data(simDat)
  
  library(snow)
  cl <- makeCluster(c("localhost","localhost"), type = "SOCK")

  path.pval <- gseaSnp(cl=cl, snp.dat=snpDat, snp.info=snp.info,
                      gene.info=gene.info, gene.set=sim.pathway, y=y,
                    snp.method="logiReg", gene.def="abs",dist=5)
  
  stopCluster(cl)
}

%\keyword{ TBD }

