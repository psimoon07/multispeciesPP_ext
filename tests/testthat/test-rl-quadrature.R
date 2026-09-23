## How the RL quadrature set is built is a data-prep decision, but it is the one most
## likely to silently corrupt the fit, so it is pinned down here.

test_that("the quadrature must match the resolution of the detection covariates", {
    ## The spec proposed one quadrature point per 100x100 m plot.  With covariates on a
    ## 25x25 m grid, that means the integral term sees plot-mean covariates while the
    ## detection margin sees cell-level ones.  The two design matrices disagree, and the
    ## fit distorts BOTH beta (to reproduce the observed covariate margin from an
    ## attenuated quadrature) and alpha (compensating).  Using the 25 m cells is exact.
    land <- sim.plot.landscape(n.plots=2500, cells.per.plot=16, within.sd=0.6, seed=3)
    coefs <- rbind(sp1=c(-2.0, 0.9, -0.6))
    lam <- sim.lambda(land, coefs)
    set.seed(808)
    w <- 3
    counts <- rpois(nrow(land), lam[,1]*w)
    RL <- list(sp1=land[rep(seq_len(nrow(land)), counts), , drop=FALSE])

    BG.cells <- land; BG.cells$quadrat.size <- w      # exact: one row per covariate cell
    BG.cent  <- sim.plot.centroids(land, w)           # one row per plot

    fit <- function(bg) multispeciesPP(~x1+x2, ~z1, PA=NULL, PO=NULL, BG=land,
                                       RL=RL, BG.RL=bg, quadrat.size.RL=bg$quadrat.size,
                                       region.size=nrow(land))
    cells <- fit(BG.cells)
    cent  <- suppressWarnings(fit(BG.cent))

    ## exact quadrature recovers the truth
    expect_lt(abs(cells$species.coef["(Intercept)",1] - coefs[1,1]), 0.1)
    expect_lt(abs(cells$species.coef["x1",1] - coefs[1,2]), 0.05)
    expect_lt(abs(cells$species.coef["x2",1] - coefs[1,3]), 0.05)

    ## the coarse quadrature inflates beta badly -- this is the real damage, and it is
    ## far worse than the intercept shift a naive Jensen argument would predict
    expect_gt(cent$species.coef["x1",1], coefs[1,2]*1.3)
    expect_gt(abs(cent$species.coef["(Intercept)",1] - coefs[1,1]), 0.25)
})

test_that("zero-detection plot-visits must stay in BG.RL", {
    ## A visit that found nothing enters the likelihood only through -integral(lambda).
    ## Filtering BG.RL down to visits that detected the species throws away every
    ## absence and inflates lambda.
    land <- sim.plot.landscape(n.plots=2500, cells.per.plot=16, within.sd=0.6, seed=3)
    coefs <- rbind(sp1=c(-4.5, 0.9, -0.6))          # rare enough that most plots are empty
    lam <- sim.lambda(land, coefs)
    set.seed(11)
    w <- 2
    counts <- rpois(nrow(land), lam[,1]*w)
    RL <- list(sp1=land[rep(seq_len(nrow(land)), counts), , drop=FALSE])
    BG.RL <- land; BG.RL$quadrat.size <- w

    per.plot <- tapply(counts, land$plot, sum)
    expect_gt(sum(per.plot == 0), 0.3*length(per.plot))   # the filter must have teeth
    detected <- land$plot %in% names(per.plot)[per.plot > 0]

    fit <- function(bg) multispeciesPP(~x1+x2, ~z1, PA=NULL, PO=NULL, BG=land,
                                       RL=RL, BG.RL=bg, quadrat.size.RL=bg$quadrat.size,
                                       region.size=nrow(land))
    complete <- fit(BG.RL)
    filtered <- suppressWarnings(fit(BG.RL[detected, ]))

    expect_lt(abs(complete$species.coef["(Intercept)",1] - coefs[1,1]), 0.15)
    ## dropping the empty visits overstates lambda by roughly e^1
    expect_gt(filtered$species.coef["(Intercept)",1] -
              complete$species.coef["(Intercept)",1], 0.5)
})
