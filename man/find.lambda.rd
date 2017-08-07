

\name{ find.lambda }
\alias{ find.lambda }
\title{ find.lambda }

\description{
  Find optimal \eqn{\lambda}{lambda} for a regularized regression 
}

\usage{
find.lambda(X, y, grp,  weights=NULL, upper=1000, intercept=TRUE, eq.eSNP=TRUE, 
            maxiter=100, tol= 0.001, method=c("AIC", "BIC"))
}

\arguments{
  \item{ X }{
    TBD
  }

  \item{ y }{
    TBD
  }

  \item{ grp }{
    TBD
  }

  \item{ weights }{
    TBD
  }

  \item{ upper }{
    TBD
  }

  \item{ intercept }{
    TBD
  }

  \item{ eq.eSNP }{
    TBD
  }

  \item{ maxiter }{
    TBD
  }

  \item{ tol }{
    TBD
  }

  \item{ method }{
    TBD
  }
}

\seealso{
  TBD
  %% \code{  \link{ help } }  
}

\details{ Searches for optimal \eqn{\lambda}{lambda} that gives an AIC value near 
  minimum. It tries a range of different \eqn{\lambda}{lambda} values when
  minimizing the penalized likelihood function, until finding a value
  that is not the maximum of the range.  If the minimum AIC value is
  on the maximum edge of the range, it retries with a larger range.  }

\value{
  return value TBD
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

