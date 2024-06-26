% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/causal_ccombat.R
\name{cb.correct.caus_cComBat}
\alias{cb.correct.caus_cComBat}
\title{Causal Conditional ComBat}
\usage{
cb.correct.caus_cComBat(
  Ys,
  Ts,
  Xs,
  match.form,
  reference = NULL,
  match.args = list(method = "nearest", exact = NULL, replace = FALSE, caliper = 0.1),
  retain.ratio = 0.05,
  apply.oos = FALSE
)
}
\arguments{
\item{Ys}{an \code{[n, d]} matrix, for the outcome variables with \code{n} samples in \code{d} dimensions.}

\item{Ts}{\code{[n]} the labels of the samples, with \code{K < n} levels, as a factor variable.}

\item{Xs}{\code{[n, r]} the \code{r} covariates/confounding variables, for each of the \code{n} samples, as a data frame with named columns.}

\item{match.form}{A formula of columns from \code{Xs}, to be passed directly to \code{\link[MatchIt]{matchit}} for subsequent matching. See \code{formula} argument from \code{\link[MatchIt]{matchit}} for details.}

\item{reference}{the name of the reference/control batch, against which to match. Defaults to \code{NULL}, which treats the reference batch as the smallest batch.}

\item{match.args}{A named list arguments for the \code{\link[MatchIt]{matchit}} function, to be used to specify specific matching strategies, where the list names are arguments and the corresponding values the value to be passed to \code{matchit}. Defaults to inexact nearest-neighbor caliper (width 0.1) matching without replacement.}

\item{retain.ratio}{If the number of samples retained is less than \code{retain.ratio*n}, throws a warning. Defaults to \code{0.05}.}

\item{apply.oos}{A boolean that indicates whether or not to apply the learned batch effect correction to non-matched samples that are still within a region of covariate support. Defaults to \code{FALSE}.}
}
\value{
a list, containing the following:
\itemize{
   \item{\code{Ys.corrected}} an \code{[m, d]} matrix, for the \code{m} retained samples in \code{d} dimensions, after correction.
   \item{\code{Ts}} \code{[m]} the labels of the \code{m} retained samples, with \code{K < n} levels.
   \item{\code{Xs}} the \code{r} covariates/confounding variables for each of the \code{m} retained samples.
   \item{\code{Model}} the fit batch effect correction model. See \code{\link[sva]{ComBat}} for details.
   \item{\code{InSample.Ids}} the ids which were used to fit the batch effect correction model.
   \item{\code{Corrected.Ids}} the ids to which batch effect correction was applied. Differs from \code{InSample.Ids} if \code{apply.oos} is \code{TRUE}.
}
}
\description{
A function for implementing the causal conditional ComBat (causal cComBat) algorithm.
This algorithm allows users to remove batch effects (in each dimension), while adjusting for known confounding
variables. It is imperative that this function is used in conjunction with domain
expertise (e.g., to ensure that the covariates are not colliders, and that the system satisfies the strong
ignorability condiiton) to derive causal conclusions. See citation for more details as to the conditions
under which conclusions derived are causal.
}
\section{Details}{

For more details see the help vignette:
\code{vignette("causal_ccombat", package = "causalBatch")}
}

\examples{
library(causalBatch)
sim <- cb.sims.sim_linear(a=-1, n=100, err=1/8, unbalancedness=3)
cb.correct.caus_cComBat(sim$Ys, sim$Ts, data.frame(Covar=sim$Xs), "Covar")

}
\references{
Eric W. Bridgeford, et al. "A Causal Perspective for Batch Effects: When is no answer better than a wrong answer?" Biorxiv (2024).

Daniel E. Ho, et al. "MatchIt: Nonparametric Preprocessing for Parametric Causal Inference" JSS (2011).

W Evan Johnson, et al. "Adjusting batch effects in microarray expression data using empirical Bayes methods" Biostatistics (2007).
}
\author{
Eric W. Bridgeford
}
