\name{multispeciesPP}
\alias{multispeciesPP}
\title{
  Fit multispecies proportional bias models.
}
\description{
  Fit the model described in the Proportional Bias Model paper, optionally combining
  up to three kinds of data source for each species: presence-absence surveys
  (\code{PA}), biased presence-only records (\code{PO}), and exhaustively censused
  plots (\code{RL}).  All sources share one species intensity
  \eqn{\lambda_k(u) = \exp(\alpha_k + \beta_k' X(u))}; only \code{PO} carries the
  observer-bias terms \eqn{\gamma_k} and \eqn{\delta}.
}
\usage{
multispeciesPP(sdm.formula, bias.formula, PA = NULL, PO = NULL, BG,
               species = union(names(PO), names(RL)),
               species.PA = species, species.PO = species,
               RL = NULL, BG.RL = NULL, species.RL = names(RL),
               quadrat.size = 1, region.size = 1, quadrat.size.RL = 1,
               start = NULL, inverse.hessian = FALSE,
               penalty.l2.sdm = 0.1, penalty.l2.bias = 0.1,
               penalty.l2.intercept = 1e-04,
               weights = rep(1, n.species * nrow(x)), control = list())
}
\arguments{
  \item{sdm.formula}{
    Formula for dependence of the species process on environmental covariates
  }
  \item{bias.formula}{
    Formula for dependence of the observer bias on covariates
  }
  \item{PA}{
    Presence-absence data, as a data frame with one row per site visit and one 0/1
    column per species.  Several presence-absence sources with different survey
    geometries may be stacked into a single frame, provided \code{quadrat.size} gives
    each row its own exposure.  \code{NULL} if there is no presence-absence data.
  }
  \item{PO}{
    Presence-only data, as a list of data frames whose names are the species.
    \code{NULL} if there is no presence-only data.
  }
  \item{BG}{
    Background data, as a data frame.  Always required: it supplies the quadrature for
    the presence-only integral and it defines the centering and scaling applied to
    every other design matrix, whether or not any \code{PO} data is present.
  }
  \item{species}{
    Names of the species.  Note that the fitted set is always
    \code{union(species.PO, species.PA, species.RL)}; a species with no data in any
    source cannot be forced into the fit by naming it here.
  }
  \item{species.PA}{
    Species contributing PA data
  }
  \item{species.PO}{
    Species contributing PO data
  }
  \item{RL}{
    Exhaustively censused point data, as a list of data frames whose names are the
    species, with one row per detected individual.  Structurally like \code{PO}, but
    with no sampling bias: no \eqn{\gamma} intercept is estimated for it and
    \code{bias.formula} does not apply, so these frames need only the covariates named
    in \code{sdm.formula}.  Must be supplied together with \code{BG.RL}.
  }
  \item{BG.RL}{
    Quadrature data frame for the RL survey extent \eqn{W_{RL}}, used to evaluate
    \eqn{\int_{W_{RL}} \lambda_k(u)\,du}.  Normally the RL plot-visits themselves,
    resolved to the covariate grid -- see the Warning section.  Must contain every
    visited location, including those where nothing was detected.
  }
  \item{species.RL}{
    Species contributing RL data.  Defaults to \code{names(RL)}.  Include a species
    only if the RL protocol would have recorded it, since the RL block assumes
    detection within a plot is exhaustive.
  }
  \item{quadrat.size}{
    Exposure of each presence-absence site visit: either a single number applied to
    every row, or a numeric vector of length \code{nrow(PA)} giving each row its own
    \eqn{A_i \cdot \Delta t_i}.  A vector may be passed unnamed (e.g. as a column of
    \code{PA}) or named by \code{rownames(PA)}; either way it is realigned internally
    after \code{model.matrix} drops rows with missing covariates.
  }
  \item{region.size}{
    Size of the study region, i.e. the measure of the domain the \code{BG} points are
    drawn from.  If background points are drawn jointly over space and time this must
    be a space-time measure, not an area.
  }
  \item{quadrat.size.RL}{
    Quadrature weight of each row of \code{BG.RL}: a single number, or a numeric
    vector of length \code{nrow(BG.RL)}.  For the usual exact scheme this is the
    space-time measure each quadrature cell represents.  To use a Monte-Carlo
    background for \eqn{W_{RL}} instead, pass random draws as \code{BG.RL} with
    \code{quadrat.size.RL = region.size.RL / nrow(BG.RL)}.
  }
  \item{start}{
    Optional argument giving initial coefficients.  These should be in
    standardized form; i.e. they should correspond to a model with the
    covariates centered and scaled.
  }
  \item{inverse.hessian}{
    Not currently used
  }
  \item{penalty.l2.sdm}{
    L2 penalty on coefficients for the species intensity
  }
  \item{penalty.l2.bias}{
    L2 penalty on coefficients for the bias
  }
  \item{penalty.l2.intercept}{
    L2 penalty on the intercepts
  }
  \item{weights}{
    Observation weights
  }
  \item{control}{
    Analogous to \code{control} for \code{glm.fit}
  }
}
\details{
  The RL block adds no coefficients.  Its rows carry \code{isPO = 0} in the species
  design and zeroes in the bias design, so its linear predictor is exactly
  \eqn{\alpha_k + \beta_k' X(u)} -- the same coefficients every other source uses.
  This is what holds \eqn{\lambda_k} common across sources, and it means RL data
  cannot move \eqn{\gamma_k} or \eqn{\delta}.  Detections enter through a compressed
  margin row and the survey extent enters through the \code{BG.RL} quadrature rows.

  Passing \code{RL = NULL} reproduces the original two-source model exactly.
}
\section{Warning}{
  \strong{Units must agree across sources.}  \eqn{\alpha_k} is shared, so
  \code{region.size}, every element of \code{quadrat.size}, and every element of
  \code{quadrat.size.RL} must be expressed in the same measure -- e.g. m^2 * day if
  \eqn{\lambda} is an intensity per square metre per day.  A mismatch does not stop
  the fit from converging; it simply returns the wrong \eqn{\alpha_k}.

  \strong{\code{BG.RL} must contain every visited location, not just the ones where
  something was found.}  A visit that detected nothing contributes to the likelihood
  only through the \eqn{-\int \lambda} term, which is carried entirely by these rows.
  Filtering \code{BG.RL} down to visits with detections discards every absence and
  inflates \eqn{\lambda_k} substantially.

  \strong{Resolve \code{BG.RL} to the covariate grid, not to the survey plot.}  One
  quadrature point per plot is exact only if \eqn{\lambda} is constant across the
  plot.  When covariates are finer than the plots -- say a 25 m raster inside a
  100 m plot -- the integral term then sees plot-mean covariates while the detection
  margin sees cell-level ones, and the fit distorts \eqn{\beta_k} as well as
  \eqn{\alpha_k}.  Use one quadrature row per covariate cell per visit, with
  \code{quadrat.size.RL} set to that cell's exposure.

  \strong{\code{species.coef["isPO", k]} is meaningful only for species in
  \code{species.PO}.}  For any other species \eqn{\gamma_k} is unidentified and the
  reported value is an artifact of the penalty, not an estimate.  This is true of the
  original two-source model as well.

  Standard errors for the intercepts are known to be incorrect.
}
\value{
  An object of class \code{multispeciesPP}, a list whose principal elements are
  \item{species.coef}{Matrix of species coefficients on the original covariate scale,
    with rows \code{(Intercept)}, the \code{sdm.formula} terms, and \code{isPO}.}
  \item{bias.coef}{Bias coefficients \eqn{\delta} on the original covariate scale.}
  \item{normalized.species.coef, normalized.bias.coef}{The same, on the standardized
    scale the model was fitted on.}
  \item{std.errs}{Standard errors; see the Warning about intercepts.}
  \item{fit.PA, fit.BG, fit.RL}{Fitted intensities at the PA sites, the background
    points, and the RL quadrature rows, each per unit of that source's exposure.}
  \item{bias.fit.BG}{Fitted sampling bias at the background points.}
  \item{species.PA, species.PO, species.RL}{The species contributing each source.}
  \item{sdm.margins.ab.RL}{Per-species column sums of the sdm design matrix over the
    RL detections, kept for diagnostics.}
}
\author{
  William Fithian
}
\seealso{
  \code{\link{predict.multispeciesPP}}, \code{\link{summary.multispeciesPP}}
}
\examples{
## Two species on a gridded landscape; lambda is constant within a cell, so every
## integral below is exact.
set.seed(1)
n.cells <- 2000
BG <- data.frame(x1 = rnorm(n.cells), x2 = rnorm(n.cells), z1 = rnorm(n.cells))
truth <- rbind(sp1 = c(-2.0, 0.8, -0.5), sp2 = c(-2.5, -0.4, 0.9))
lambda <- exp(cbind(1, BG$x1, BG$x2) \%*\% t(truth))
region.size <- n.cells          # one unit of exposure per cell

## Presence-only: thinned by exp(gamma_k + delta * z1)
gamma <- c(sp1 = 1.2, sp2 = 1.5); delta <- 0.6
bias <- exp(outer(BG$z1 * delta, gamma, "+"))
PO <- lapply(1:2, function(k)
    BG[rep(1:n.cells, rpois(n.cells, lambda[, k] * bias[, k])), ])
names(PO) <- rownames(truth)

## Two presence-absence sources with different transect geometries, stacked
sites <- sample(n.cells, 1200)
PA <- BG[sites, ]
PA$quadrat.size <- rep(c(0.5, 8), length.out = length(sites))   # A_i * dt_i
for (k in 1:2)
    PA[[rownames(truth)[k]]] <- rbinom(length(sites), 1,
                                       1 - exp(-lambda[sites, k] * PA$quadrat.size))

## Exhaustively censused plots: detections, plus a quadrature row for EVERY visit
plots <- sample(n.cells, 800)
BG.RL <- BG[plots, ]
BG.RL$quadrat.size <- 20
RL <- lapply(1:2, function(k)
    BG[plots[rep(seq_along(plots), rpois(length(plots), lambda[plots, k] * 20))], ])
names(RL) <- rownames(truth)

fit <- multispeciesPP(~ x1 + x2, ~ z1, PA = PA, PO = PO, BG = BG,
                      RL = RL, BG.RL = BG.RL,
                      quadrat.size = PA$quadrat.size,
                      quadrat.size.RL = BG.RL$quadrat.size,
                      region.size = region.size)
fit$species.coef     # compare with t(truth) and gamma
fit$bias.coef        # compare with delta
}
\keyword{models}
\keyword{spatial}
