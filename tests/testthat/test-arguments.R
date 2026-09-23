## Argument handling for the new optional sources.

small <- function() {
    land <- sim.landscape(n.cells=1200, seed=2)
    lam <- sim.lambda(land, rbind(sp1=c(-2, .6, -.4)))
    set.seed(6)
    rl <- sim.RL(land, lam, plots=sample(nrow(land), 600), weight=15)
    list(land=land, lam=lam, RL=rl$RL, BG.RL=rl$BG.RL)
}

test_that("RL and BG.RL must be supplied together", {
    d <- small()
    expect_error(multispeciesPP(~x1+x2, ~z1, BG=d$land, RL=d$RL),
                 "both be supplied, or both be NULL")
    expect_error(multispeciesPP(~x1+x2, ~z1, BG=d$land, BG.RL=d$BG.RL),
                 "both be supplied, or both be NULL")
})

test_that("a call with no species data is rejected", {
    d <- small()
    expect_error(multispeciesPP(~x1+x2, ~z1, BG=d$land), "no species have data")
})

test_that("a badly sized quadrat.size.RL is rejected", {
    d <- small()
    expect_error(multispeciesPP(~x1+x2, ~z1, BG=d$land, RL=d$RL, BG.RL=d$BG.RL,
                                quadrat.size.RL=c(1,2,3)),
                 "scalar or a numeric vector of length nrow\\(BG.RL\\)")
})

test_that("a species whose only source is RL but which was never detected warns", {
    d <- small()
    expect_warning(multispeciesPP(~x1+x2, ~z1, BG=d$land,
                                  RL=list(sp1=d$land[0, , drop=FALSE]), BG.RL=d$BG.RL,
                                  quadrat.size.RL=d$BG.RL$quadrat.size,
                                  region.size=nrow(d$land)),
                   "no RL detections")
})

test_that("species is inferred from PO and RL names together", {
    ## A species whose only source is RL must not be silently dropped from the fit.
    land <- sim.landscape(n.cells=1500, seed=8)
    truth <- rbind(sp1=c(-2,.6,-.4), sp2=c(-2.5,-.3,.7))
    lam <- sim.lambda(land, truth)
    set.seed(21)
    PO <- sim.PO(land, lam, gamma=c(sp1=1, sp2=1.2), delta=.4)
    rl <- sim.RL(land, lam, plots=sample(nrow(land), 700), weight=15)

    ## sp2 has RL only; sp1 has PO only.  Both must survive into the fit.
    fit <- suppressMessages(
        multispeciesPP(~x1+x2, ~z1, PA=NULL, PO=PO["sp1"], species.PO="sp1", BG=land,
                       RL=rl$RL["sp2"], BG.RL=rl$BG.RL, species.RL="sp2",
                       quadrat.size.RL=rl$BG.RL$quadrat.size, region.size=nrow(land)))

    expect_equal(sort(colnames(fit$species.coef)), c("sp1","sp2"))
    expect_equal(fit$species.RL, "sp2")
    expect_equal(fit$species.PO, "sp1")

    ## and sp2's intensity is recovered from its RL data alone
    expect_lt(abs(fit$species.coef["(Intercept)","sp2"] - truth["sp2",1]), 0.2)
    expect_lt(abs(fit$species.coef["x1","sp2"] - truth["sp2",2]), 0.1)
})

test_that("the units reminder fires only when it is relevant", {
    d <- small()
    expect_message(multispeciesPP(~x1+x2, ~z1, BG=d$land, RL=d$RL, BG.RL=d$BG.RL,
                                  quadrat.size.RL=d$BG.RL$quadrat.size,
                                  region.size=nrow(d$land)),
                   "same space-time units")

    land <- sim.landscape(n.cells=1200, seed=2)
    lam <- sim.lambda(land, rbind(sp1=c(-2,.6,-.4)))
    set.seed(6)
    PO <- sim.PO(land, lam, gamma=c(sp1=1), delta=.4)
    PA <- sim.PA(land, lam, sites=sample(nrow(land), 400), quadrat.size=2)
    expect_no_message(multispeciesPP(~x1+x2, ~z1, PA=PA, PO=PO, BG=land,
                                     quadrat.size=2, region.size=nrow(land)))
})
