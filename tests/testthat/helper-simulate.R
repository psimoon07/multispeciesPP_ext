## Simulation helpers shared by the test suite.
##
## Space is discretized into `n.cells` cells of equal area `cell.area`, each carrying
## its own covariates.  lambda_k is constant within a cell, so every integral below is
## exact and the tests compare against true coefficients rather than a quadrature
## approximation of them.  All exposures (cell.area, quadrat sizes, RL plot weights,
## region.size) are in the same units, as the model requires.

sim.landscape <- function(n.cells=4000, cell.area=1, seed=1) {
    set.seed(seed)
    data.frame(x1=rnorm(n.cells), x2=rnorm(n.cells), z1=rnorm(n.cells),
               cell.area=cell.area, row.names=paste0("c",seq_len(n.cells)))
}

## A landscape whose cells are grouped into survey plots, as the RL design is: a
## 100x100 m plot holds 16 cells of the 25x25 m covariate grid.  Covariates are
## spatially autocorrelated -- a plot-level mean plus a smaller within-plot
## deviation -- which is what makes the choice of quadrature resolution a question
## of degree rather than a pathology.  `within.sd` controls how much covariate
## variation a single plot-centroid quadrature point would throw away.
sim.plot.landscape <- function(n.plots=3000, cells.per.plot=16, cell.area=1,
                               within.sd=0.6, seed=1) {
    set.seed(seed)
    n.cells <- n.plots*cells.per.plot
    plot.id <- rep(seq_len(n.plots), each=cells.per.plot)
    between.sd <- sqrt(max(0, 1 - within.sd^2))   # keep total variance at 1
    jitter <- function() rnorm(n.plots, sd=between.sd)[plot.id] + rnorm(n.cells, sd=within.sd)
    data.frame(x1=jitter(), x2=jitter(), z1=jitter(),
               cell.area=cell.area, plot=plot.id,
               row.names=paste0("c",seq_len(n.cells)))
}

## Collapse a plot landscape to one centroid row per plot -- i.e. spec 2b applied
## literally, one quadrature point per 100x100 m plot.
sim.plot.centroids <- function(land, weight.per.cell) {
    cen <- aggregate(cbind(x1,x2,z1) ~ plot, data=land, FUN=mean)
    cen$quadrat.size <- weight.per.cell * as.vector(table(land$plot)[as.character(cen$plot)])
    rownames(cen) <- paste0("p",cen$plot)
    cen
}

## lambda_k per unit area at every cell, for a matrix of coefficients whose
## rows are species and whose columns are (intercept, x1, x2).
sim.lambda <- function(land, coefs) {
    X <- cbind(1, land$x1, land$x2)
    lam <- exp(X %*% t(coefs))
    dimnames(lam) <- list(rownames(land), rownames(coefs))
    lam
}

## Presence-only: thin the true process by b_k(u) = exp(gamma_k + delta * z1(u)),
## then return one row per detected point.
sim.PO <- function(land, lam, gamma, delta) {
    b <- exp(outer(land$z1 * delta, gamma, "+"))
    counts <- matrix(rpois(length(lam), lam * b * land$cell.area), nrow(lam))
    setNames(lapply(seq_len(ncol(lam)), function(k)
        land[rep(seq_len(nrow(land)), counts[,k]), , drop=FALSE]),
        colnames(lam))
}

## Presence-absence: `sites` indexes cells, `quadrat.size` is per-site exposure.
## Detection is exhaustive within the quadrat, hence cloglog.
sim.PA <- function(land, lam, sites, quadrat.size) {
    q <- rep(quadrat.size, length.out=length(sites))
    out <- land[sites, , drop=FALSE]
    for(k in colnames(lam))
        out[[k]] <- rbinom(length(sites), 1, 1 - exp(-lam[sites,k] * q))
    out$quadrat.size <- q
    rownames(out) <- paste0("pa", seq_len(nrow(out)))
    out
}

## Rote-Liste: exhaustive census of `plots` (cell indices) with exposure `weight`.
## Returns the per-species detection frames and the quadrature frame, which contains
## EVERY visited plot -- including those where nothing was found, since a
## non-detection enters the likelihood only through the -integral(lambda) term.
sim.RL <- function(land, lam, plots, weight) {
    w <- rep(weight, length.out=length(plots))
    counts <- matrix(rpois(length(plots)*ncol(lam), lam[plots,] * w), nrow=length(plots))
    RL <- setNames(lapply(seq_len(ncol(lam)), function(k)
        land[plots[rep(seq_along(plots), counts[,k])], , drop=FALSE]),
        colnames(lam))
    BG.RL <- land[plots, , drop=FALSE]
    BG.RL$quadrat.size <- w
    rownames(BG.RL) <- paste0("rl", seq_along(plots))
    list(RL=RL, BG.RL=BG.RL)
}
