

\name{ grass }
\alias{ grass }
\title{ grass }

\description{
    This function calculates the p-value for the association between a pathway and
    disease risk using the GRASS algorithm. GRASS summarizes the genetic structure
    by singular value decomposition for each gene as eigenSNPs and uses a novel form
    of regularized regression technique, termed group ridge regression, to select 
    representative eigenSNPs for each gene and assess their joint association with 
    disease risk. 
}

\usage{
    grass(cl=NULL, snp.dat, snp.info, gene.info, gene.set, y, weights=NULL, 
      gene.def=c("abs", "rel"), dist=20000,  k=1, B=100, nominal.p=FALSE)
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
    inverse probability that the subject is selected.  
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
  
  \item{ nominal.p }{A logical value indicating if p-values should be calculated based on normal
    approximation of null statistics. The default is FALSE, which means p-values are calculated
    by permutation. When computation load is heavy, one can use a small number of 
    permutations \link{B} and set {nominal.p = TRUE}. 
  }


}

\seealso{
    \code{  \link{ aligator } } , \code{  \link{ gseaSnp } }, \code{  \link{ plinkSet } }   
}

\details{
    The GRASS algorithm first summarizes the genetic structure (e.g., LD) for each gene as
    independent eigenSNPs via principal component analysis. Next, GRASS uses regularized 
    regression technique with a novel form of penalty function, group ridge penalty, to 
    estimate the association of eigenSNPs with disease risk. The group ridge regularized 
    regression applies a lasso penalty within a gene to select representative eigenSNPs 
    while simultanously using a ridge penalty among the genes to shrink the gene estimate. 
    In this way, efficient estimation is achieved even with high-dimensional data and the 
    pathway statistics are not dominated by a single gene or a single SNP. A pathway
    statistic is calculated from the penalized log odds ratios of eigenSNPs in the
    pathway. By permuting disease status, recomputing the null pathway statistics and comparing
    them with the observed ones, we obtain the p-value of disease association for each pathway.
}

\value{
  Returns a numerical vector of p-values, corresponding to the pathways in \link{gene.set}.  
}

\references{
  Chen LS, Hutter CM, Potter JD, Liu Y, Prentice RL, Peters U, Hsu L (2010). Insights 
  into colon cancer etiology using a regularized approach to gene set analysis of GWAS
  data.  American Journal of Human Genetics, 86(6): 860-871.
}

\author{
  Lin S. Chen <lche2@fhcrc.org>
}

%\note{
%  TBD
%}

\examples{
  library(SNPath)
  data(simDat)
  
  library(snow)
  cl <- makeCluster(c("localhost","localhost"), type = "SOCK")
  clusterEvalQ(cl, {
    library(SNPath)
  })

  path.pval <- grass(cl=cl, snp.dat=snpDat, snp.info=snp.info, gene.info=gene.info,
                   gene.set=sim.pathway, y=y, gene.def="rel", B=1000)
  
  stopCluster(cl)
}

%\keyword{ TBD }

