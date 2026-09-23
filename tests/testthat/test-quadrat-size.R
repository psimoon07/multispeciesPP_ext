## Change 1: quadrat.size as a per-row exposure.  ST and ZK have different transect
## geometries, so stacking them into one PA frame under a single scalar offset is wrong.

two.source.PA <- function() {
    land <- sim.landscape(n.cells=4000, seed=7)
    coefs <- rbind(sp1=c(-2.0, 0.8, -0.5), sp2=c(-2.5, -0.4, 0.9))
    lam <- sim.lambda(land, coefs)
    set.seed(4242)
    st <- sim.PA(land, lam, sites=sample(nrow(land), 900), quadrat.size=0.5)  # small transects
    zk <- sim.PA(land, lam, sites=sample(nrow(land), 900), quadrat.size=8)    # large transects
    PA <- rbind(st, zk)
    rownames(PA) <- paste0("pa", seq_len(nrow(PA)))
    list(land=land, coefs=coefs, PA=PA,
         PO=sim.PO(land, lam, gamma=c(sp1=1.2, sp2=1.5), delta=0.6))
}

test_that("a per-row quadrat.size recovers alpha where a single scalar does not", {
    d <- two.source.PA()
    fit <- function(q) multispeciesPP(~x1+x2, ~z1, PA=d$PA, PO=d$PO, BG=d$land,
                                      quadrat.size=q, region.size=nrow(d$land))
    per.row <- fit(d$PA$quadrat.size)
    scalar  <- fit(mean(d$PA$quadrat.size))

    truth <- d$coefs[,1]
    err.per.row <- abs(per.row$species.coef["(Intercept)",] - truth)
    err.scalar  <- abs(scalar$species.coef["(Intercept)",]  - truth)

    expect_lt(max(err.per.row), 0.25)       # recovers ground truth
    expect_true(all(err.scalar > err.per.row))
    expect_gt(min(err.scalar), 0.3)         # and the scalar fit is visibly biased
})

test_that("quadrat.size accepts a bare column or a vector named by rownames(PA)", {
    d <- two.source.PA()
    fit <- function(q) multispeciesPP(~x1+x2, ~z1, PA=d$PA, PO=d$PO, BG=d$land,
                                      quadrat.size=q, region.size=nrow(d$land))
    unnamed <- fit(d$PA$quadrat.size)
    named   <- fit(setNames(d$PA$quadrat.size, rownames(d$PA)))
    expect_equal(unnamed$all.coef, named$all.coef)

    ## a scalar must still behave exactly like the equivalent constant vector
    const.scalar <- fit(3)
    const.vector <- fit(rep(3, nrow(d$PA)))
    expect_equal(const.scalar$all.coef, const.vector$all.coef)
})

test_that("quadrat.size stays aligned when model.matrix drops NA covariate rows", {
    ## PA.good.rows is a character vector of rownames and can be a strict subset of
    ## rownames(PA).  Subsetting quadrat.size by name rather than position is the bug
    ## this test exists to catch.
    d <- two.source.PA()
    PA.na <- d$PA
    PA.na$x1[c(3, 17, 400, 1200)] <- NA
    PA.pre <- PA.na[!is.na(PA.na$x1), ]
    expect_lt(nrow(PA.pre), nrow(PA.na))

    fit <- function(P) multispeciesPP(~x1+x2, ~z1, PA=P, PO=d$PO, BG=d$land,
                                      quadrat.size=P$quadrat.size, region.size=nrow(d$land))
    expect_equal(fit(PA.na)$all.coef, fit(PA.pre)$all.coef)
})

test_that("a badly sized quadrat.size is rejected rather than silently recycled", {
    d <- two.source.PA()
    expect_error(multispeciesPP(~x1+x2, ~z1, PA=d$PA, PO=d$PO, BG=d$land,
                                quadrat.size=c(1,2), region.size=nrow(d$land)),
                 "scalar or a numeric vector of length nrow\\(PA\\)")
})
