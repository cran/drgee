\encoding{latin1}
\name{drgee}
\alias{drgee}
\alias{coef.drgee}
\alias{vcov.drgee}
\alias{print.drgee}
\alias{summary.drgee}
\alias{print.summary.drgee}
\title{Doubly Robust Generalized Estimating Equations}
\description{
  \code{drgee} performs outcome nuisance model based estimation,
  exposure nuisance model based estimation or doubly robust estimation
  given symbolic representations of an outcome nuisance model and
  an exposure nuisance model.
}
\usage{
drgee(oformula, eformula, iaformula = formula(~1),
      olink = c("identity", "log", "logit"),
      elink = c("identity", "log", "logit"),
      estimationMethod = c("dr", "obe", "ebe"),
      data = NULL, rootFinder = findRoots,
      clusterid = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{oformula}{
    An expression or formula for the outcome nuisance model.
    The outcome is identified as the response in
    this formula. Therefore, a LHS is required for all choices for
    \code{estimationMethod}.
  }
  \item{eformula}{
    An expression or formula for the exposure nuisance model. The
    exposure is identified as the response in this formula. Therefore, a
    LHS is required for all choices for \code{estimationMethod}.
}
\item{iaformula}{
  An expression or formula where the RHS should contain the variables
  that "interact" (i.e. are supposed to be multiplied with) with the
  exposure in the main model. "1" will always added. Default value is no
  interactions,
  i.e. \code{iaformula=formula(~1)}.
}
\item{olink}{
  A character string naming the link function in the outcome nuisance
  model. Has to be \code{"identity"}, \code{"log"} or
  \code{"logit"}. Default is \code{"identity"}.
}
\item{elink}{
  A character string naming the link function in the exposure nuisance
  model. Has to be \code{"identity"}, \code{"log"} or
  \code{"logit"}. Default is \code{"identity"}.
}
\item{estimationMethod}{
  A character string naming the desired estimation method. Choose
  \code{"obe"} for outcome nuisance model based estimation,
  \code{"ebe"} for exposure nuisance model based estimation or
  \code{"dr"} for doubly robust estimation. Default is \code{"dr"}.
}
\item{data}{
   A data frame or environment containing the variables appearing in \code{iaformula},
    \code{oformula} and \code{eformula}. Default is \code{NULL} in which
    case data are expected to be found in the environment of the
    \code{oformula} argument.
}
\item{rootFinder}{
  A function to solve a system of non linear equations. Default
  is \code{findRoots}.
}
\item{clusterid}{
  A optional character string naming a cluster-defining variable in the
  \code{data} argument.
}
\item{\dots}{
  Further arguments to be passed to the functions \code{geeFit} and \code{rootFinder}.
}
}
\details{
  \code{drgee} estimates the parameter \eqn{\beta}{beta} in a main
  model \eqn{g\{E(Y|A,L)\}=\beta^T \{A\cdot X(L)\}+Q(L)}{g\{E(Y|A,L)\}=beta^T (A * X(L))+Q(L)\}},
  where \eqn{L} is a vector of nuisance variables and
  \eqn{X(L)} and \eqn{Q(L)} are functions of \eqn{L}. Note that \eqn{A
    \cdot X(L)}{A*X(L)} should be interpreted as a columnwise
  multiplication and that \eqn{X(L)} will always contain a column of 1's.
  Given a specification of an outcome nuisance model \eqn{Q(L)=\gamma^T
    V(L)}{Q(L)=gamma^T V(L)} (where \eqn{V(L)} is a function of \eqn{L})
  outcome nuisance model based estimation can be
  performed. Alternatively, leaving \eqn{Q(L)} unspecified and
  using an exposure nuisance model \eqn{h\{E(A|L)\}=\alpha^T
    Z(L)}{h\{E(Y|L)\}=alpha^T Z(L)} (where \eqn{h} is a link
  function and \eqn{Z(L)} is a function of \eqn{L}) exposure
  nuisance model based estimation can be performed.  When \eqn{g} is
  \eqn{logit}, the exposure nuisance model is required be
  of the form \eqn{logit\{E(A|Y=0,L)\}=\alpha^T Z(L)}{logit\{E(A|Y=0,L)\}=alpha^T Z(L)}.
  In this case the exposure needs to binary. Given both
  an outcome and an exposure nuisance model, doubly robust
  estimation can be performed. Doubly robust estimation
  gives a consistent estimate of the parameter \eqn{\beta}{beta} when
  either the outcome nuisance model or the exposure nuisance model
  is correctly specified, not necessarily both.




  Usage is best explained through an example. Suppose that we are
  interested in the parameter vector \eqn{\beta_0}{beta_0} and
  \eqn{\beta_1}{beta_1} in a main model \eqn{logit\{E(Y|A,L_1,L_2)\}=\beta_0 A + \beta_1
    A \cdot L_1 + Q(L_1,L_2)}{logit\{E(Y|A,L_1,L_2)\}=beta_0 A + beta_1 A *
    L_1 +
    Q(L_1,L_2)} where \eqn{L_1} and \eqn{L_2} are nuisance variables and
  \eqn{Q(L_1,L_2)} is some (unspecified) function of \eqn{L_1} and \eqn{L_2}.

  To adjust for  \eqn{L_1} and \eqn{L_2}, we can use an outcome
  nuisance model \eqn{Q(L_1,L_2;\gamma)=\gamma_0 + \gamma_1
  L_1}{Q(L_1,L_2;gamma)=gamma_0 + gamma_1 L_1} or an
  exposure nuisance model \eqn{logit\{E(A|Y=0,L_1,L_2)\}=\alpha_0+\alpha_1
  L1+\alpha_2 L2}{logit\{E(A|Y=0,L_1,L_2)\}=alpha_0+alpha_1
  L1+alpha_2 L2} to calculate estimates of \eqn{\beta_0}{beta_0} and \eqn{\beta_1}{beta_1}
  in the main model.

  We specify the outcome nuisance model as \code{oformula=Y~L_1}
  and \code{olink="logit"}. The exposure nuisance model is specified as
  \code{eformula=A~L_1+L_2} and \code{elink="logit"}.
  Since the outcome \eqn{Y} and the exposure \eqn{A} are
  identified as the LHS of \code{oformula} and \code{eformla}
  respectively and since the outcome link is specified in the
  \code{olink} argument,
  the only thing left to specify for the main model is the
  (multiplicative) interactions \eqn{X(L)=(1,L_1)^T}. This is done as
  \code{iaformula=~L_1}, since \eqn{1} is always included in \eqn{X(L)}.
  We can then perform outcome or exposure nuisance model based
  estimation or doubly robust estimation by setting
  \code{estimationMethod} to \code{"obe"}, \code{"ebe"} or \code{"dr"} respectively.

  When \code{estimationMethod="obe"}, the RHS of \code{eformula} will be ignored
  with a warning message.

  When \code{estimationMethod="ebe"}, the RHS of \code{oformula} will be ignored with a warning message.

  Outcome nuisance model based estimation is implemented for generalized
  estimating equation models with the identity, log or logit link and
  independent observations. The estimated coefficients are identical to
  those obtained with \code{glm}, but since no
  distributional assumptions are made, robust variance is calculated.

  When exposure nuisance model based estimation or doubly robust estimation
  estimation is chosen with \code{olink="logit"} the exposure link will be
  changed to \code{"logit"} with a warning message.

  Robust variance for the estimated parameter is calculated
  using \code{robVcov}. A cluster robust variance is calculated when
  a character string naming a cluster variable is
  supplied in the \code{clusterid} argument.

  \code{drgee} calls \code{geeData} to create a \code{geeData} object
  containing the elements needed in the calculations. The estimation of
  the coefficents in the main model is performed by \code{obeFit},
  \code{ebeFit} or \code{drFit}.

  For exposure nuisance model based estimation when \eqn{g} is the identity
  or log link, see Robins et al. (1992).

  For doubly robust estimation when \eqn{g} is the identity
  or log link, see Robins (1999). For doubly robust estimation when
  \eqn{g} is the logit link, see Tchetgen et al. (2010).

  This package was inspired by the STATA package \code{drglm} described
  in Orsini et al. (2013). Basically, it also provides the same functionality.

}
\value{
  \code{drgee} return an object of class \code{drgee} containing:

  \item{coefficients }{Estimates of the parameters in the main model.}
  \item{vcov }{Robust variance of the parameter estimates.}
  \item{optim.object }{An estimation object returned from the function specified
  	in the \code{rootFinder}, if this function is called.}
  \item{call }{The matched call.}
  \item{geeData }{The \code{geeData} object used in the calculations.}
  \item{estimationMethod }{The value of the input argument \code{estimationMethod}.}

  The class methods \code{coef} and \code{vcov} can be used to extract
  the estimated parameters and their covariance matrix from a
  \code{drgee} object. \code{summary.drgee} produces a summary of the
  calculations.
}
\seealso{
  \code{\link{obeFit}} for outcome nuisance model based estimation,
  \code{\link{ebeFit}} for exposure nuisance model based
  estimation, \code{\link{drFit}} for doubly robust estimation
  of the parameters in the main model, \code{\link{drgeeData}} for data
  preparation and \code{\link{findRoots}} for nonlinear equation
  solving, \code{\link{robVcov}} for estimation of variance.
}
\author{
Johan Zetterqvist, Arvid \enc{Sjölander}{Sjolander}
}

\examples{

## Doubly robust estimation when
## the main model is
## E(Y|A,L1,L2)-E(Y|A=0,L1,L2)=beta0*A+beta1*A*L1
## and the outcome nuisance model is
## E(Y|A=0,L1,L2)=gamma0+gamma1*L1+gamma2*L2
## and the exposure nuisance model is
## E(A|Y=0,L1,L2)=expit(alpha0+alpha1*L1+alpha2*l2)

library(drgee)

expit<-function(x) exp(x)/(1+exp(x))

n<-5000

# nuisance
l1<-rnorm(n, mean = 0, sd = 1)
l2<-rnorm(n, mean = 0, sd = 1)

beta0<-1.5
beta1<-1
gamma0<--1
gamma1<--2
gamma2<-2
alpha0<-1
alpha1<-5
alpha2<-3

# Exposure
a<-rbinom(n,1,expit(alpha0 + alpha1*l1 + alpha2*l2))
# Outcome
y<-rnorm(n,beta0*a + beta1*a*l1 + gamma0 + gamma1*l1 + gamma2*l2,sd=1)

data<-data.frame(y,a,l1,l2)

## outcome nuisance model misspecified and
## exposure nuisance model correctly specified

# Doubly robust estimation
dr.est <- drgee(y~l1,a~l1+l2,~l1,"identity","logit","dr",data)
summary(dr.est)

# Outcome nuisance model based estimation
obe.est <- drgee(y~l1,a~1,~l1,"identity","logit","obe",data)
summary(obe.est)

# Exposure based estimation
ebe.est <- drgee(y~1,a~l1+l2,~l1,"identity","logit","ebe",data)
summary(ebe.est)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
\keyword{robust}

\references{
  Orsini N., Belocco R., \enc{Sjölander}{Sjolander} A. (2013), Doubly
  Robust Estimation in Generalized Linear Models, \emph{Stata Journal},
  \bold{13}, 1, pp.185-205

  Robins J.M., Mark S.D., Newey W.K. (1992), Estimating Exposure
  Effects by Modelling the Expectation of Exposure Conditional
  on Confounders, \emph{Biometrics}, \bold{48}, 479--495

  Robins JM (1999), Robust Estimation in Sequentially Ignorable
  Missing Data and Causal Inference Models, \emph{Proceedings of the
  American Statistical Association Section on Bayesian Statistical
  Science}, pp. 6--10

  Tchetgen E.J.T., Robins J.M., Rotnitzky A. (2010), On Doubly Robust
  Estimation in a Semiparametric Odds Ratio Model, \emph{Biometrika}, \bold{97},
  1, 171--180
 }
