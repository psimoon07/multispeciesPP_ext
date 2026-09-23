# Extending `multispeciesPP` to Four Spatio-Temporal Data Sources

## Context

`multispeciesPP` (Fithian, 2014) fits a joint multispecies IPP model from exactly
two data sources per species:

- `PA`: a single presence-absence data frame, fit as `binomial(link="cloglog")`
  rows, with one **global scalar** `quadrat.size` offset for every row.
- `PO`: a list of per-species presence-only data frames, fit via a Poisson
  margin/quadrature trick against one shared background set `BG`, with a
  shared bias slope `δ` and species-specific bias intercept `γ_k`.

We now have **four** sources per species (see the project's derivation
document, §6): PO, Stichproben (ST), Zufallskoordinaten (ZK), and Rote-Liste
plots (RL). ST and ZK are structurally identical to the existing `PA` type
(Bernoulli/cloglog, no bias term). RL is a *new* likelihood type: an unbiased,
exhaustive point pattern (own coordinate per detection, no `γ`/`δ`), restricted
to its own survey extent — i.e. "PO with `b_k(u) ≡ 1`".

This document specifies the two code changes needed, in order of size. All
`§` references below are to `R/multispeciesPP.R` unless noted.

---

## Change 1 (small): `quadrat.size` must be a per-row exposure, not one scalar

**Where:** `multispeciesPP.R`, the `offset` construction loop (currently
around the `for(k in 1:n.species)` block, line `offk[1 + n.species + (1:n.PA)]
<- log(quadrat.size)`).

**Why:** ST and ZK have different transect geometries (`A_ST = 200 m²` vs.
`A_ZK ≈ 452 m²`), each surveyed for a fixed duration `Δt` (one day). The
correct offset per row is `log(q_i)` where `q_i = A_i · Δt_i`, not a single
constant for the whole `PA` data frame.

**Fix:**

1. Accept `quadrat.size` as either a scalar (backward-compatible) **or** a
   numeric vector of length `nrow(PA)`.
2. **Critical bug to avoid:** `PA.good.rows` (and thus `n.PA`) can be a strict
   subset of `nrow(PA)` — rows get dropped by `model.matrix()` if any
   covariate is `NA`. If `quadrat.size` is passed in as a vector aligned to
   the *original* `PA` row order, it must be subset/reordered by
   `PA.good.rows` exactly like `sdm.PA.model.matrix` is, **before** being used
   in the offset assignment. Concretely: expect the caller to supply
   `quadrat.size` as a column of `PA` itself (e.g. `PA$quadrat.size`), or as a
   named vector keyed by `rownames(PA)`, and inside the function do
   `quadrat.size <- quadrat.size[PA.good.rows]` right after `PA.good.rows` is
   computed.
3. No other change is needed — `offk[...] <- log(quadrat.size)` already
   broadcasts correctly over a vector once it's the right length and order.

**Usage after the fix:** stack Stichproben and Zufallskoordinaten site-visits
into one `PA` data frame (species columns as 0/1, as today), plus a
`quadrat.size` column giving each row's `A_i · Δt_i` in the same units. No
other part of the source-blending logic changes.

**Test:** with two synthetic PA sources of different known quadrat sizes,
confirm recovered `α_k, β_k` match ground truth (they should differ visibly
from a fit that ignores the per-row exposure and uses a single scalar).

---

## Change 2 (real addition): a third likelihood block for RL data

**Why RL doesn't fit the existing two types:** unlike ST/ZK, RL records exact
coordinates for every detected individual within exhaustively-searched
100×100 m plots — so it's a genuine point-process likelihood, not
Bernoulli. Unlike PO, it has **no sampling bias**: `γ_RL` doesn't exist and
`δ` doesn't apply. Mathematically (see derivation §6.3, Eq. RL):

```
ℓ_k,RL = Σ_i log λ_k(u_RL,i) − ∫_{W_RL} λ_k(u) du
```

which is structurally identical to the PO block's margin/quadrature
machinery, **except**: (a) no bias covariates enter (`z ≡ 0` for every RL
row), (b) no `γ_RL` intercept is estimated, and (c) the quadrature/background
set must be restricted to `W_RL` (the union of the 189 RL plot-visits), not
the full domain `W` used for PO.

### 2a. New function arguments

Add to `multispeciesPP()`:

```r
RL = NULL,            # list of per-species data frames, like PO
BG.RL = NULL,          # background/quadrature data frame for the RL domain
region.size.RL = 1,    # |W_RL|; ignored if using exact quadrature (2b)
species.RL = species   # which species have RL data
```

### 2b. Recommended quadrature: use the RL visits themselves

Rather than drawing Monte-Carlo background points for `BG.RL` (as PO does
against `BG`), **use the actual RL plot-visit list as the quadrature set**,
each with weight `A_RL · Δt`. This is exact (Berman & Turner, 1992) because
RL plots are exhaustively censused — there is no need to approximate `∫_{W_RL}
λ_k du` via random background draws when the true survey design already
tiles `W_RL` at known, fixed locations/dates. This sidesteps needing a
separate `region.size.RL` / `n.BG.RL` bookkeeping step entirely: every RL
plot-visit becomes one quadrature row with offset `log(A_RL · Δt)`,
regardless of whether that visit detected the species or not.

(If a genuinely separate background set is preferred instead, follow the
existing `PO`/`BG` pattern with `region.size.RL` and `offset = log(region.size.RL)
- log(n.BG.RL)`, restricted to spatio-temporal draws inside `W_RL`.)

### 2c. Matrix construction changes

Mirror the existing PO block (`sdm.margins.ab`, `sdm.PO.model.matrices`,
`PO.good.rows`, the `x`/`z`/`y`/`offset` append loop, `which.BG`-style index
vector) but as a **separate, parallel block**:

1. Compute `sdm.margins.ab.RL` (column sums of the sdm design matrix over
   RL detections, per species) exactly as `sdm.margins.ab` is computed for
   PO — **do not** add these into `sdm.margins.ab` itself; keep them as
   their own margin rows.
2. Do **not** compute or accumulate anything into `bias.margins` from RL
   data. The RL rows' `z`-columns are all zero.
3. Do **not** add an `"isRL"` intercept column analogous to PO's `"isPO"`
   column in `sdm.margins.ab`/`sdm.coef`. RL's linear predictor is exactly
   `α_k + β_k^T X(u)`, reusing the *same* `α_k, β_k` columns as PO/PA/ST/ZK
   — this is what keeps `λ_k` identical across all four sources, per the
   derivation's joint-likelihood requirement (§6.4).
4. Append the RL quadrature rows (§2b) as `poisson()`-family rows, with `z =
   0` throughout, exactly analogous to how BG rows are appended for PO
   (same `row.families` pattern, new block index).
5. Update `row.families` to include the new block type, and update the
   `n.species*(...)` bookkeeping in `abcd.from.all.coef` /
   `all.coef.from.abcd` **only if** a new coefficient block was added — with
   the design above, no new coefficients are introduced, so these two
   functions should not need to change.
6. Update the score-check consistency block (`score.check.gamma`,
   `score.check.bias`) if RL margins participate in the gradient check —
   since RL doesn't touch `γ` or `δ`, its only contribution to any score
   check would be through `α_k, β_k`, which are already checked implicitly
   via the existing PA/PO score equations; no new score-check variable
   should be needed, but confirm convergence diagnostics still pass with RL
   rows included.

### 2d. `predict.multispeciesPP.R`

No change needed: RL introduces no new coefficients, so `object$species.coef`
and `object$bias.coef` retain their existing meaning and shape.

**Test:** simulate RL-only data under a known `λ_k` with **no** bias applied
(direct realization of the true IPP, no thinning), fit with only the RL block
active (PO/PA/ST/ZK empty), and confirm recovered `α_k, β_k` are unbiased —
in particular confirm the fit does **not** try to estimate a `γ` or `δ`
effect from RL rows (both should stay exactly at their PA/PO-block-driven
values, unaffected by RL data).

---

## Change 3 (no code — data prep only): time-varying covariates

The package has no notion of time; every row's covariates are just whatever
columns are in `PA`/`PO`/`RL`/`BG`. To reuse it for the spatio-temporal model:

- Every PA/ST/ZK/RL visit and every PO record needs its covariate columns
  (14-day trailing temperature/precipitation, `sin`/`cos` seasonal terms,
  centered year) computed **at that record's actual date**, not a static
  site value.
- Background/quadrature points (`BG` for PO; `BG.RL` if used instead of
  §2b's exact-quadrature approach) must be drawn jointly over **space and
  time** — i.e. each background point needs an independently-sampled date
  from the relevant observation window, matching what the manuscript's
  Methods §2.3 already does for the 2M PO background points.
- No package-level change is required for this; it is purely a matter of
  how the `PA`, `PO`, `RL`, and `BG` data frames are constructed before
  being passed in.

---

## Summary of deliverables for this repo

1. `R/multispeciesPP.R`: vectorize + correctly align `quadrat.size` (Change 1).
2. `R/multispeciesPP.R`: add the `RL`/`BG.RL`/`region.size.RL` arguments and
   the new margin+quadrature block, with `z ≡ 0` and no new intercept
   (Change 2).
3. No changes needed to `predict.multispeciesPP.R`, `summary.multispeciesPP.R`,
   `block.glm.fit.R`, or `block.projection.R`.
4. Add/update `man/multispeciesPP.Rd` to document the new arguments.
5. Add unit tests per the "Test" notes above (quadrat-size recovery test;
   RL-only unbiasedness test; combined four-source fit sanity check against
   the existing two-source fit as a special case with `RL = NULL`).
