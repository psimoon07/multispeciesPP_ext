## Change 2: the RL likelihood block.  RL is an exhaustive, unbiased point pattern --
## "PO with b_k(u) == 1" -- so it shares alpha/beta with every other source but must
## touch neither gamma nor delta.

rl.setup <- function(seed=4242, n.cells=4000) {
    land <- sim.landscape(n.cells=n.cells, seed=7)
    coefs <- rbind(sp1=c(-2.0, 0.8, -0.5), sp2=c(-2.5, -0.4, 0.9))
    lam <- sim.lambda(land, coefs)
    set.seed(seed)
    rl <- sim.RL(land, lam, plots=sample(nrow(land), 1500), weight=20)
    list(land=land, coefs=coefs, lam=lam, RL=rl$RL, BG.RL=rl$BG.RL,
         region.size=nrow(land)*land$cell.area[1])
}

test_that("an RL-only fit recovers alpha and beta without inventing a bias effect", {
    d <- rl.setup()
    fit <- multispeciesPP(~x1+x2, ~z1, PA=NULL, PO=NULL, BG=d$land,
                          RL=d$RL, BG.RL=d$BG.RL,
                          quadrat.size.RL=d$BG.RL$quadrat.size,
                          region.size=d$region.size)

    expect_lt(max(abs(fit$species.coef["(Intercept)",] - d$coefs[,1])), 0.15)
    expect_lt(max(abs(fit$species.coef["x1",] - d$coefs[,2])), 0.1)
    expect_lt(max(abs(fit$species.coef["x2",] - d$coefs[,3])), 0.1)

    ## RL rows carry isPO = 0 in x and z = 0, so neither gamma nor delta sees them at
    ## all: both sit exactly where the ridge penalty alone puts them.
    expect_equal(unname(fit$species.coef["isPO",]), c(0, 0))
    expect_equal(unname(fit$bias.coef), 0)
})

test_that("RL data for one species leaves the other species untouched", {
    ## Regression test for the NA mask on RL quadrature rows.  Without it, a species
    ## with no RL data still absorbs integral(lambda) over W_RL and its alpha sinks.
    d <- rl.setup()
    lam <- d$lam
    set.seed(99)
    PO <- sim.PO(d$land, lam, gamma=c(sp1=1.2, sp2=1.5), delta=0.6)
    PA <- sim.PA(d$land, lam, sites=sample(nrow(d$land), 900), quadrat.size=3)

    base <- function(...) multispeciesPP(~x1+x2, ~z1, PA=PA, PO=PO, BG=d$land,
                                         quadrat.size=PA$quadrat.size,
                                         region.size=d$region.size, ...)
    ## Both fits halve a step on the way in and warn "stopped at boundary value".  That
    ## is pre-existing IRLS behaviour -- the unmodified function warns identically on
    ## this data and returns identical coefficients -- and the score checks inside
    ## multispeciesPP() are what guarantee the fit still lands on the optimum.
    without <- suppressWarnings(base())
    with.sp1 <- suppressWarnings(
        base(RL=d$RL["sp1"], BG.RL=d$BG.RL, species.RL="sp1",
             quadrat.size.RL=d$BG.RL$quadrat.size))

    expect_equal(with.sp1$species.coef[, "sp2"], without$species.coef[, "sp2"],
                 tolerance=1e-3)
    expect_false(isTRUE(all.equal(with.sp1$species.coef[, "sp1"],
                                  without$species.coef[, "sp1"], tolerance=1e-3)))
})

test_that("a four-source joint fit recovers every coefficient", {
    d <- rl.setup()
    set.seed(99)
    PO <- sim.PO(d$land, d$lam, gamma=c(sp1=1.2, sp2=1.5), delta=0.6)
    st <- sim.PA(d$land, d$lam, sites=sample(nrow(d$land), 900), quadrat.size=0.5)
    zk <- sim.PA(d$land, d$lam, sites=sample(nrow(d$land), 900), quadrat.size=8)
    PA <- rbind(st, zk); rownames(PA) <- paste0("pa", seq_len(nrow(PA)))

    fit <- multispeciesPP(~x1+x2, ~z1, PA=PA, PO=PO, BG=d$land,
                          RL=d$RL, BG.RL=d$BG.RL,
                          quadrat.size=PA$quadrat.size,
                          quadrat.size.RL=d$BG.RL$quadrat.size,
                          region.size=d$region.size)

    expect_lt(max(abs(fit$species.coef["(Intercept)",] - d$coefs[,1])), 0.15)
    expect_lt(max(abs(fit$species.coef["x1",] - d$coefs[,2])), 0.1)
    expect_lt(max(abs(fit$species.coef["x2",] - d$coefs[,3])), 0.1)
    expect_lt(max(abs(fit$species.coef["isPO",] - c(1.2, 1.5))), 0.2)
    expect_lt(abs(fit$bias.coef - 0.6), 0.1)

    ## RL adds no coefficients, so the returned shapes are unchanged
    expect_equal(rownames(fit$species.coef), c("(Intercept)","x1","x2","isPO"))
    expect_equal(length(fit$bias.coef), 1L)
    expect_equal(dim(fit$fit.RL), c(nrow(d$BG.RL), 2L))
})

test_that("RL detections need only the sdm covariates, not the bias covariates", {
    d <- rl.setup()
    RL.nobias <- lapply(d$RL, function(df) df[, c("x1","x2","cell.area"), drop=FALSE])
    BG.RL.nobias <- d$BG.RL[, c("x1","x2","quadrat.size"), drop=FALSE]
    expect_silent({
        fit <- suppressMessages(
            multispeciesPP(~x1+x2, ~z1, PA=NULL, PO=NULL, BG=d$land,
                           RL=RL.nobias, BG.RL=BG.RL.nobias,
                           quadrat.size.RL=BG.RL.nobias$quadrat.size,
                           region.size=d$region.size))
    })
    expect_lt(max(abs(fit$species.coef["(Intercept)",] - d$coefs[,1])), 0.15)
})
