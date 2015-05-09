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
drgeeData(outcome, exposure,
          oformula, eformula, iaformula = formula(~1),
          olink = c("identity", "log", "logit"),
          elink = c("identity", "log", "logit"),
          data, estimation.method = c("dr", "o", "e"),
          cond = FALSE, clusterid)
}

\arguments{
  \item{outcome }{
    The outcome as a variable or as a character string naming a variable in the
    \code{data} argument. If it is not found in the \code{data}
    argument, it will be searched for in the calling frame. If missing,
    the outcome is assumed to be the response of \code{oformula}.
  }

  \item{exposure }{
    The exposure as a variable or as a character string naming a variable in the
    \code{data} argument. If it is not found in the \code{data}
    argument, it will be searched for in the calling frame. If missing,
    the outcome is assumed to be the response of \code{eformula}.
  }

  \item{oformula }{
    An expression or formula for the outcome nuisance model. The outcome
    is identified as the response in this formula.
  }

  \item{eformula }{
    An expression or formula for the exposure nuisance model. The
    exposure is identified as the response in this formula.
  }

  \item{iaformula }{
    An expression or formula where the RHS should contain the variables
    that "interact" (i.e. are supposed to be multiplied with) with the
    exposure in the main model to create the terms associated with the
    parameters of interest. "1" will always added. Default value is no
    interactions, i.e. \code{formula(~1)}.
  }

  \item{olink }{
    A character string naming the link function in the outcome nuisance
    model. Have to be \code{"identity"}, \code{"log"} or
    \code{"logit"}. Default is \code{"identity"}.
  }

  \item{elink }{
    A character string naming the link function in the exposure nuisance
    model. Have to be \code{"identity"}, \code{"log"} or
    \code{"logit"}. Default is \code{"identity"}. When
    \code{olink="logit"} this is replaced by \code{"logit"}.
  }

  \item{data}{
    A data frame or environment containing the variables in \code{iaformula},
    \code{oformula} and \code{eformula}. If missing, variables are
    expected to be found in the calling frame.
  }

  \item{estimation.method}{
    A character string naming the desired estimation method. Choose
    \code{"o"} for O-estimation,
    \code{"e"} for E-estimation or
    \code{"dr"} for DR-estimation. Default is \code{"dr"}.
  }

  \item{cond }{
    A logical value indicating whether the nuisance models should have
    cluster-specific intercepts. If \code{cond=TRUE} the design matrices
    for the nuisance models do not have an intercept. Requires a
    \code{clusterid} argument.
  }

  \item{clusterid}{
    A cluster-defining variable or a character string naming a cluster-defining variable in the
    \code{data} argument. If it is not found in the \code{data}
    argument, it will be searched for in the calling frame. If missing,
    each observation will be considered to be a separate cluster. This
    argument is required when \code{cond = TRUE}.
  }

}
\details{
  \code{drgeeData} is called by \code{drgee} and \code{gee} to extract
  data from a \code{data.frame} or \code{environment}
  object. The data can then be used to for O-estimation, E-estimation or
  DR-estimation. \code{drgeeData} uses
  \code{model.frame} and \code{model.matrix} to remove incomplete
  observations and to convert factors to dummy variables. It also
  performs check the supplied data for errors or inconsistencies.

  The class method \code{summary.drgeeData} produces strings for the
  formulas with terms referring to the columns in the produced design
  matrices.
}
\value{
  \code{drgee.data } returns an object of class \code{drgeeData}
  containing
  \item{y }{The outcome matrix.}
  \item{outname }{A string for the name of the outcome.}
  \item{a }{The exposure matrix.}
  \item{expname }{A string for the name of the exposure.}
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
  \item{clustname }{A string for the name of the cluster defining variable.}
  \item{olink }{A character string naming the link function in the outcome nuisance
  model.}
  \item{elink }{A character string naming the link function in the outcome nuisance
    model.}
  \item{cond }{A logical value indicating whether cluster-specific intercepts should
  be assumed. If \code{TRUE}, the is no column for the intercept in
  \code{v} and \code{z}. Outcome concordant will also be removed.}
  \item{oterms }{The \code{terms} object corresponding to the outcome
  nuisance model.}
  \item{eterms }{The \code{terms} object corresponding to the exposure
  nuisance model.}

  All matrix elements have rows corresponding to the complete observations
  in the original data.
}
\author{
  Johan Zetterqvist, Arvid \enc{Sjölander}{Sjolander}
}
\seealso{
\code{\link{drgee}}, \code{\link{gee}}, \code{\link{model.frame}} and \code{\link{model.matrix}}.
}
\keyword{models}
\keyword{regression}
\keyword{cluster}
