% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/causal_ccombat.R
\name{cb.correct.aipw_cComBat}
\alias{cb.correct.aipw_cComBat}
\title{Augmented Inverse Probability Weighting Conditional ComBat}
\usage{
cb.correct.aipw_cComBat(
  Ys,
  Ts,
  Xs,
  aipw.form,
  covar.out.form = NULL,
  retain.ratio = 0.05
)
}
\arguments{
\item{Ys}{an \code{[n, d]} matrix, for the outcome variables with \code{n} samples in \code{d} dimensions.}

\item{Ts}{\code{[n]} the labels of the samples, with at most two unique batches.}

\item{Xs}{\code{[n, r]} the \code{r} covariates/confounding variables, for each of the \code{n} samples, as a data frame with named columns.}

\item{aipw.form}{A covariate model, given as a formula. Applies for the estimation of propensities for the AIPW step.}

\item{covar.out.form}{A covariate model, given as a formula. Applies for the outcome regression step of the \code{ComBat} algorithm. Defaults to \code{NULL}, which re-uses \code{aipw.form} for the covariate/outcome model.}

\item{retain.ratio}{If the number of samples retained is less than \code{retain.ratio*n}, throws a warning. Defaults to \code{0.05}.}
}
\value{
a list, containing the following:
\itemize{
   \item{\code{Ys.corrected}} an \code{[m, d]} matrix, for the \code{m} retained samples in \code{d} dimensions, after correction.
   \item{\code{Ts}} \code{[m]} the labels of the \code{m} retained samples, with \code{K < n} levels.
   \item{\code{Xs}} the \code{r} covariates/confounding variables for each of the \code{m} retained samples.
   \item{\code{Model}} the fit batch effect correction model.
   \item{\code{Corrected.Ids}} the ids to which batch effect correction was applied.
}
}
\description{
A function for implementing the AIPW conditional ComBat (AIPW cComBat) algorithm.
This algorithm allows users to remove batch effects (in each dimension), while adjusting for known confounding
variables. It is imperative that this function is used in conjunction with domain
expertise (e.g., to ensure that the covariates are not colliders, and that the system could be argued to satisfy the
ignorability condition) to derive causal conclusions. See citation for more details as to the conditions
under which conclusions derived are causal.
}
\details{
Note: This function is experimental, and has not been tested on real data. It has only been tested with simulated data with binary (0 or 1) exposures.
}
\section{Details}{

For more details see the help vignette:
\code{vignette("causal_ccombat", package = "causalBatch")}
}

\examples{
library(causalBatch)
sim <- cb.sims.sim_linear(a=-1, n=100, err=1/8, unbalancedness=2)
cb.correct.aipw_cComBat(sim$Ys, sim$Ts, data.frame(Covar=sim$Xs), "Covar")

}
\references{
Eric W. Bridgeford, et al. "A Causal Perspective for Batch Effects: When is no answer better than a wrong answer?" Imaging Neuroscience (2025).

W Evan Johnson, et al. "Adjusting batch effects in microarray expression data using empirical Bayes methods" Biostatistics (2007).
}
\author{
Eric W. Bridgeford
}
