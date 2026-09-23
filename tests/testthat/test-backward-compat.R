## The RL block inserts a new row block into x/z and shifts n.sites.  This test is the
## guard that none of that index surgery changed the original two-source model.

test_that("a two-source fit is bit-for-bit identical to the pre-extension function", {
    land <- sim.landscape(n.cells=3000, seed=11)
    coefs <- rbind(sp1=c(-2.0, 0.8, -0.5), sp2=c(-2.5, -0.4, 0.9))
    lam <- sim.lambda(land, coefs)
    set.seed(99)
    PO <- sim.PO(land, lam, gamma=c(sp1=1.2, sp2=1.5), delta=0.6)
    PA <- sim.PA(land, lam, sites=sample(nrow(land), 800), quadrat.size=3)
    region.size <- nrow(land)*land$cell.area[1]

    old <- multispeciesPP.baseline(~x1+x2, ~z1, PA=PA, PO=PO, BG=land,
                                   quadrat.size=3, region.size=region.size)
    new <- multispeciesPP(~x1+x2, ~z1, PA=PA, PO=PO, BG=land,
                          quadrat.size=3, region.size=region.size)

    expect_equal(new$all.coef, old$all.coef)
    expect_equal(new$normalized.all.coef, old$normalized.all.coef)
    expect_equal(new$std.errs, old$std.errs)
    expect_equal(new$species.coef, old$species.coef)
    expect_equal(new$bias.coef, old$bias.coef)
    expect_equal(new$fit.PA, old$fit.PA)
    expect_equal(new$fit.BG, old$fit.BG)
})

test_that("RL=NULL leaves every pre-existing return value in place", {
    land <- sim.landscape(n.cells=1500, seed=5)
    lam <- sim.lambda(land, rbind(sp1=c(-2, .6, -.4)))
    set.seed(3)
    PO <- sim.PO(land, lam, gamma=c(sp1=1), delta=.5)
    PA <- sim.PA(land, lam, sites=sample(nrow(land), 400), quadrat.size=2)
    fit <- multispeciesPP(~x1+x2, ~z1, PA=PA, PO=PO, BG=land,
                          quadrat.size=2, region.size=nrow(land))

    expect_true(all(c("species.coef","bias.coef","fit.PA","fit.BG","bias.fit.BG",
                      "linear.fit.PA","linear.fit.BG","all.coef","std.errs")
                    %in% names(fit)))
    ## predict() keys its default species list off colnames(fit.PA)
    expect_equal(colnames(fit$fit.PA), "sp1")
    expect_equal(nrow(predict(fit, land[1:10,])), 10L)
})
