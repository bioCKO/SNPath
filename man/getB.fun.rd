

\name{ getB.fun }
\alias{ getB.fun }
\title{ getB.fun }

\description{
  TBD
}

\usage{
getB.fun(gseaLis, y0=NULL, weights0=NULL, lambda=NULL, nmsiz=100, intercept=TRUE, 
         eq.eSNP=TRUE,tol=0.001,maxiter=100, method="AIC", 
         return.lambda=TRUE)
}

\arguments{
  \item{ gseaLis }{
    %% index into global gsea.obj list object
  }

  \item{ y0 }{
    If non-null, use y0 for disease status.
  }

  \item{ weights0 }{
  }

  \item{ lambda }{
    If null, compute lambda using find.lambda
  }

  \item{ nmsiz }{
    maximum number of genes used at a time
  }

  \item{ intercept }{
    passed to find.lambda, and logistic.ridge.c
  }

  \item{ eq.eSNP }{
     passed to find.lambda, and logistic.ridge.c
  }

  \item{ tol }{
    passed to find.lambda, and logistic.ridge.c
  }

  \item{ maxiter }{
    passed to find.lambda, and logistic.ridge.c
  }

  \item{ method }{
    passed to find.lambda
  }

  \item{ return.lambda }{
    if true, returns list containing beta vector and lambda vector
  }
}

\seealso{
  TBD
  %% \code{  \link{ help } }  
}

\details{
  If needed, split the \eqn{X}{X} matrix of characteristic vectors 
}

\value{
  Returns \eqn{\beta}{beta} list or list containing \eqn{\beta}{beta} list and \eqn{\lambda}{lambda} value if return.lambda is true.
}

\references{
  TBD
}

%\author{
%  TBD
%}

\note{
  TBD
}

\examples{
  #TBD
}

\keyword{ TBD }

