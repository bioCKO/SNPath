


\name{ calc.fun }
\alias{ calc.fun }
\title{ calc.fun }

\description{
    This function calculates the marginal association p-values and test statistics
    for individual SNPs.
}

\usage{
    calc.fun(cl=NULL, snp.dat, y, weights=NULL, snp.method=c("logiReg","chiSq") )
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

}

\seealso{
  \code{  \link{ aligator } } , \code{  \link{ gseaSnp } }, \code{  \link{ grass } }, 
  \code{ \link{ plinkSet }}
}


\value{
  Returns a list of two vectors: 
  \item{pval}{A numerical vector of p-values for each SNP.}
  \item{stat}{A numerical vector of test statistics for each SNP.}  
}




\examples{
  data(simDat)
  
  library(snow)
  cl <- makeCluster(c("localhost","localhost"), type = "SOCK")

  snp.pval <- calc.fun(cl=cl, snp.dat=snpDat, y=y, weights=NULL, snp.method="logiReg" )

  ## weighted analysis
  ## require the library Zelig and survey
  
  library(Zelig)
  library(survey)
  weights <- runif(length(y), 0.5, 1.5)  ## simulate some weights
  snp.pval.w <- calc.fun(cl=cl, snp.dat=snpDat, y=y, weights=weights, snp.method="logiReg" )

}


