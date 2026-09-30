## predict() must rebuild the design from the terms stored at fit time.  Terms such as
## ns() or poly() fix their knots / centres on BG; re-deriving them from newdata gives
## a different basis, i.e. silently wrong predictions -- or an error when newdata has
## too few distinct values to place knots at all.

spline.fit <- function() {
    land <- sim.landscape(n.cells=3000, seed=21)
    coefs <- rbind(sp1=c(-2.0, 0.8, -0.5), sp2=c(-2.5, -0.4, 0.9))
    lam <- sim.lambda(land, coefs)
    set.seed(7)
    PO <- sim.PO(land, lam, gamma=c(sp1=1.2, sp2=1.5), delta=0.6)
    PA <- sim.PA(land, lam, sites=sample(nrow(land), 800), quadrat.size=3)
    fit <- multispeciesPP(~x1+splines::ns(x2, df=3), ~z1, PA=PA, PO=PO, BG=land,
                          quadrat.size=3, region.size=nrow(land)*land$cell.area[1])
    list(fit=fit, land=land)
}

## The sdm part of the fit's own BG linear predictor: what predict(sdm=TRUE) must return.
sdm.part <- function(fit, rows) {
    fit$linear.fit.BG[rows,,drop=FALSE] - fit$linear.bias.fit.BG[rows,,drop=FALSE] -
        rep(fit$normalized.species.coef["isPO",], each=length(rows))
}

test_that("predict on a subset of BG reproduces the fitted linear predictor", {
    d <- spline.fit()
    rows <- 1:40
    expect_equal(predict(d$fit, d$land[rows,]), sdm.part(d$fit, rows),
                 ignore_attr=TRUE)
})

test_that("predict works when newdata cannot place the knots itself", {
    d <- spline.fit()
    one <- d$land[c(5,5,5),]            # a single distinct value of x2
    rownames(one) <- paste0("r", 1:3)
    p <- predict(d$fit, one)
    expect_equal(dim(p), c(3L, 2L))
    expect_equal(unname(p[1,]), unname(sdm.part(d$fit, 5)[1,]))
})

test_that("a fit without stored terms still predicts from its formula", {
    d <- spline.fit()
    old <- d$fit
    old$sdm.terms <- old$bias.terms <- NULL
    ## on the full BG the formula re-derives the same knots, so the two agree
    expect_equal(predict(old, d$land), predict(d$fit, d$land))
})
