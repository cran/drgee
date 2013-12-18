\encoding{latin1}
\name{drgeeData}
\alias{drgeeData}
\alias{summary.drgeeData}
\alias{print.summary.drgeeData}
\alias{drgeeData}

\title{
  Extracting Variables and Model Matrices for Generalized Estimating equations
}
\description{
  Given a main model, an outcome nuisance model and an
  exposure nuisance model \code{drgeeData} extracts the model
  variables and matrices from a \code{data.frame} or an \code{environment}
  object. It also performs some data cleaning and error checking.
}
\usage{
drgeeData(oformula, eformula, iaformula = formula(~1),
              olink = c("identity", "log", "logit"),
              elink = c("identity", "log", "logit"),
              data = NULL, clusterid = NULL)
}

\arguments{
  \item{oformula}{
    An expression or formula for the outcome nuisance model. The outcome
  is identified as the response in this formula.
  }
  \item{eformula}{
    An expression or formula for the exposure nuisance model. The
    exposure is identified as the response in this formula.
  }
  \item{iaformula}{
    An expression or formula where the RHS should contain the variables
    that "interact" (i.e. are supposed to be multiplied with) with the
    exposure in the main model to create the terms associated with the
    parameters of interest. "1" will always added. Default value is no
    interactions, i.e. \code{formula(~1)}.
  }
  \item{olink}{
    A character string naming the link function in the outcome nuisance
  model. Have to be \code{"identity"}, \code{"log"} or
  \code{"logit"}. Default is \code{"identity"}.
  }
  \item{elink}{
     A character string naming the link function in the exposure nuisance
  model. Have to be \code{"identity"}, \code{"log"} or
  \code{"logit"}. Default is \code{"identity"}. When
    \code{olink="logit"} this is replaced by \code{"logit"}.
  }
  \item{data}{
    A data frame or environment containing the variables in \code{iaformula},
    \code{oformula} and \code{eformula}. Default is \code{NULL} in which
    case the variables will be expected to be found in the environment of the
    \code{oformula} argument.
  }
  \item{clusterid}{
      A optional character string naming a cluster-defining variable in the
  \code{data} argument.
  }
}
\details{
  \code{drgeeData} is called by \code{drgee} and uses the
  \code{oformula}, \code{eformula} and \code{iaformula} arguments to
  extract data from a \code{data.frame} or \code{environment}
  object. The data can then be used to for outcome/exposure nuisance model
  based estimation or doubly robust estimation by calling \code{obeFit},
  \code{ebeFit} or \code{drFit} respectively. \code{drgeeData} uses
  \code{model.frame} and \code{model.matrix} to remove incomplete
  observations and to convert factors to dummy variables. Factor outcome
  is only allowed when \code{olink="logit"} and then only when the
  factor has two levels. There are no restrictions on the exposure (RHS of
  \code{eformula}) or on factors appearing on the RHS of \code{oformula},
  \code{eformula} and \code{iaformula}.

  The class method \code{summary.drgeeData} produces strings for the
  formulas with terms referring to the columns in the produced design
  matrices.
}
\value{
  \code{drgeeData } returns an object of class \code{drgeeData}
  containing
  \item{y }{The outcome matrix.}
  \item{a }{The exposure matrix.}
  \item{x }{The matrix of of interactions defined in \code{iaformula}.
  	   This matrix will always contain a column with 1's.}
  \item{ax }{The matrix of elementwise product(s) of \code{a} and
  each column in \code{x}.}
  \item{v }{The matrix of terms in the outcome nuisance model.}
  \item{z }{The matrix of terms in the exposure nuisance model.}
  \item{yx }{The matrix of elementwise product(s) of \code{y} and each
  column in \code{x}. }
  \item{id }{A factor defining clusters. For independent observations,
  the number of levels equals the number of complete observations.}
  \item{olink }{ A character string naming the link function in the outcome nuisance
  model.}
  \item{elink }{ A character string naming the link function in the outcome nuisance
  model.}

  All matrix elements have rows corresponding to the complete observations
  in the original data.
}
\author{
  Johan Zetterqvist, Arvid \enc{Sjölander}{Sjolander}
}
\seealso{
\code{\link{drgee}}, \code{\link{model.frame}} and \code{\link{model.matrix}}.
}
\keyword{models}
