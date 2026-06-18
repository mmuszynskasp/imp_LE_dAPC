###############################################################################
## Project:  Why Life Expectancy Improvements Slowed (dAPC / LE)
## Aim:      LE improvements, Monte Carlo confidence bands for the LE-improvement series AND the period change-decomposition.
##           Change-decomposition over fixed periods (LE_change_decomp.csv): average annual change in LE split into drift / period / cohort contribution, with CIs.
## Inputs:   <models.dir>/fits_bs10w.rds   (total, from 02)
##           <models.dir>/fits_indiv10w.rds(baseline, from 04)
##           <data.dir>/lexis.csv          (to map model rows to A, P, Y)
## Outputs:  <models.dir>/LE_series_ci.csv
##           <models.dir>/LE_change_decomp.csv  (+ _wide.csv) shap_drift, shap_period, shap_cohort (Total; rolling Shapley, sum to dAPC)
## Depends:  dplyr, tidyr, MASS, Epi (for the fits)
###############################################################################
rm(list = ls())

library(dplyr)
library(tidyr)
library(MASS)
library(Epi)
library(xtable)

data.dir   <- #"your/input/directory"  
models.dir <- #"your/output/directory"  


B    <- 2000          # number of draws
TOL  <- 1e-2          # relative tolerance for the drift validation

fits_tot  <- readRDS(file.path(models.dir, "fits_bs10w.rds"))
fits_base <- readRDS(file.path(models.dir, "fits_indiv10w.rds"))
lexis     <- read.table(file.path(data.dir, "lexis.csv"), sep = ",", header = TRUE)

## periods for the change-decomposition (Total rates). The last is shorter (4y)
## -> its CI is wider; that is expected, not instability.
periods <- data.frame(
  period = c("1983-1993", "1993-2003", "2003-2013", "2013-2017"),
  t0     = c(1983, 1993, 2003, 2013),
  t1     = c(1993, 2003, 2013, 2017),
  stringsAsFactors = FALSE
)

## period lengths (years) used to length-weight the 1983-2017 summary row
wlen <- c("1983-1993" = 10, "1993-2003" = 10, "2003-2013" = 10, "2013-2017" = 4)

## ---- helpers ---------------------------------------------------------------

## Partial LE over ascending ages, ax = 0.5, open last interval. M = nAge x nCol.
pLE_matrix <- function(M) {
  K <- nrow(M)
  qx <- M / (1 + 0.5 * M); qx[K, ] <- 1
  px <- 1 - qx
  lx <- matrix(1, K, ncol(M))
  if (K > 1) for (k in 2:K) lx[k, ] <- lx[k - 1, ] * px[k - 1, ]
  dx <- matrix(0, K, ncol(M))
  if (K > 1) dx[1:(K - 1), ] <- lx[1:(K - 1), ] - lx[2:K, ]
  dx[K, ] <- lx[K, ]
  Lx <- matrix(0, K, ncol(M))
  if (K > 1) Lx[1:(K - 1), ] <- lx[2:K, ] + 0.5 * dx[1:(K - 1), ]
  Lx[K, ] <- lx[K, ] / M[K, ]
  colSums(Lx)
}

## Y-weighted slope extractor for the period block: logdrift = dvec %*% b_period.
make_drift_vec <- function(Xp, u, w) {
  Sw <- sum(w); Swu <- sum(w * u); Swuu <- sum(w * u * u)
  denom <- Sw * Swuu - Swu^2
  (Sw * colSums((w * u) * Xp) - Swu * colSums(w * Xp)) / denom
}

## Same window/conventions as pLE_matrix, plus restricted e-dagger and entropy.
## ed = remaining temp. LE at age of death (ax=0.5); H = edag/LE. Returns 3 x M
## with rows LE, edag, H. LE row equals pLE_matrix(M) exactly (cross-check below).
pLE_edag_matrix <- function(M) {
  K <- nrow(M); Mc <- ncol(M)
  qx <- M / (1 + 0.5 * M); qx[K, ] <- 1
  px <- 1 - qx
  lx <- matrix(1, K, Mc)
  if (K > 1) for (k in 2:K) lx[k, ] <- lx[k - 1, ] * px[k - 1, ]
  dx <- matrix(0, K, Mc)
  if (K > 1) dx[1:(K - 1), ] <- lx[1:(K - 1), ] - lx[2:K, ]
  dx[K, ] <- lx[K, ]
  Lx <- matrix(0, K, Mc)
  if (K > 1) Lx[1:(K - 1), ] <- lx[2:K, ] + 0.5 * dx[1:(K - 1), ]
  Lx[K, ] <- lx[K, ] / M[K, ]
  Tx <- apply(Lx, 2, function(z) rev(cumsum(rev(z))))      # K x M
  ex <- Tx / lx                                            # remaining temp LE at exact age
  ed <- 0.5 * rbind(ex[-1, , drop = FALSE], 0) + 0.5 * ex  # at age of death
  LE <- colSums(Lx); edag <- colSums(dx * ed)
  rbind(LE = LE, edag = edag, H = edag / LE)
}

## 10-year rolling average of annual changes; returns complete rows only.
impr_mat <- function(LE, years) {
  d1 <- LE[-1, , drop = FALSE] - LE[-nrow(LE), , drop = FALSE]
  dy <- years[-1]; T1 <- nrow(d1)
  cs <- rbind(0, apply(d1, 2, cumsum))
  idx <- 10:T1
  out <- (cs[idx + 1, , drop = FALSE] - cs[idx - 9, , drop = FALSE]) / 10
  list(impr = out, years = dy[idx])
}

## tidy a (point, draws) pair into estimate + 95% band
band_df <- function(years, point_vec, draw_mat, series, rates, cc, sx) {
  qs <- apply(draw_mat, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  tibble(country = cc, sex = sx, year = as.numeric(years),
         series = series, rates = rates,
         estimate = as.numeric(point_vec), lower = qs[1, ], upper = qs[2, ])
}

## within-period average annual change; x is a (year x B) matrix OR a year vector
period_change <- function(x, years, t0, t1) {
  i0 <- match(t0, years); i1 <- match(t1, years)
  if (is.matrix(x)) (x[i1, ] - x[i0, ]) / (t1 - t0)
  else              (x[i1]   - x[i0])   / (t1 - t0)
}

## ---- core: one fit -> LE trajectories (draws + point) for the 4 series ------
process_fit <- function(m, dat, B, seed) {
  b <- coef(m$Model); nm <- names(b)
  ## quasi-Poisson dispersion correction: the Poisson vcov is too small under
  ## over-dispersion (phi up to ~11 for the large populations), so scale it by
  ## phi = sum(Pearson^2)/df before drawing. max(.,1) never narrows below Poisson.
  phi <- sum(residuals(m$Model, type = "pearson")^2) / m$Model$df.residual
  phi <- max(phi, 1)
  V   <- phi * vcov(m$Model)
  ia <- which(startsWith(nm, "MA"))
  ip <- which(startsWith(nm, "MPr"))
  ic <- which(startsWith(nm, "MCr"))
  X  <- model.matrix(m$Model)
  stopifnot(nrow(X) == nrow(dat), length(ia) + length(ip) + length(ic) == length(b))
  Xa <- X[, ia, drop = FALSE]; Xp <- X[, ip, drop = FALSE]; Xc <- X[, ic, drop = FALSE]
  
  refP <- as.numeric(m$Ref["Per"])
  u <- dat$P - refP; w <- dat$Y
  dvec <- make_drift_vec(Xp, u, w)
  
  ## validate the drift extractor against apc.fit's reported drift
  drift_computed <- exp(sum(dvec * b[ip]))
  drift_reported <- as.numeric(m$Drift[1, 1])
  if (abs(drift_computed - drift_reported) / drift_reported > TOL)
    stop(sprintf("drift mismatch: computed %.5f vs reported %.5f", drift_computed, drift_reported))
  
  ageint <- floor(dat$A); yearint <- floor(dat$P)
  yrs <- sort(unique(yearint)); ags <- sort(unique(ageint))
  cellf <- factor(paste(yearint, ageint, sep = "_"))
  denom <- rowsum(w, cellf)[, 1]
  
  le_from_lr <- function(LR) {                 # LR: n x M log-rates -> year x M LE
    rate  <- exp(LR)
    numer <- rowsum(rate * w, cellf)
    cellrate <- numer / denom
    LE <- matrix(NA_real_, length(yrs), ncol(LR))
    for (iy in seq_along(yrs)) {
      rn <- paste(yrs[iy], ags, sep = "_")
      LE[iy, ] <- pLE_matrix(cellrate[rn, , drop = FALSE])
    }
    LE
  }
  
  ## draws
  set.seed(seed)
  Bd <- MASS::mvrnorm(B, b, V)
  Ea <- Xa %*% t(Bd[, ia, drop = FALSE])
  Ep <- Xp %*% t(Bd[, ip, drop = FALSE])
  Ec <- Xc %*% t(Bd[, ic, drop = FALSE])
  ld <- as.vector(Bd[, ip, drop = FALSE] %*% dvec)   # logdrift per draw
  U  <- outer(u, ld)
  
  draws <- list(
    dAPC = le_from_lr(Ea + Ep + Ec),
    dAP  = le_from_lr(Ea + Ep),
    dAC  = le_from_lr(Ea + U + Ec),
    dA   = le_from_lr(Ea + U)
  )
  
  ## MLE point (drift anchored at refP)
  ea <- Xa %*% b[ia]; ep <- Xp %*% b[ip]; ec <- Xc %*% b[ic]
  um <- u * sum(dvec * b[ip])
  point <- list(
    dAPC = le_from_lr(matrix(ea + ep + ec, ncol = 1))[, 1],
    dAP  = le_from_lr(matrix(ea + ep,      ncol = 1))[, 1],
    dAC  = le_from_lr(matrix(ea + um + ec, ncol = 1))[, 1],
    dA   = le_from_lr(matrix(ea + um,      ncol = 1))[, 1]
  )
  ## restricted e-dagger / entropy for the drift-only (dA) point series (H1)
  edagH_from_lr <- function(LR) {
    cellrate <- rowsum(exp(LR) * w, cellf) / denom
    t(sapply(seq_along(yrs), function(iy)
      pLE_edag_matrix(cellrate[paste(yrs[iy], ags, sep = "_"), , drop = FALSE])))
  }
  dA_edagH <- edagH_from_lr(matrix(ea + um, ncol = 1))   # years x 3: LE, edag, H
  stopifnot(isTRUE(all.equal(dA_edagH[, 1], point$dA)))  # LE must match pLE_matrix
  
  list(years = yrs, draws = draws, point = point,
       drift = c(computed = drift_computed, reported = drift_reported),
       edagH = dA_edagH)
}

## ---- main loop over country x sex ------------------------------------------
keys   <- intersect(names(fits_tot), names(fits_base))
out    <- list()    # improvement-series bands
decomp <- list()    # period change-decomposition (sub-periods)
decomp_draws <- list()   # per-draw component vectors (for the summary row CI)
decomp_pts   <- list()   # per-period POINT components (for the weighted point)
edagH_out <- list()   # drift-only e-dagger and entropy (H1 absolute-vs-relative)

for (key in keys) {
  parts <- strsplit(key, "_")[[1]]; cc <- parts[1]; sx <- parts[2]
  dat <- lexis %>% filter(country == cc, sex == sx) %>%
    transmute(A = Age, P = Year, D = round(Dx), Y = round(Exp)) %>% as.data.frame()
  
  rt <- process_fit(fits_tot[[key]],  dat, B, seed = 1000 + match(key, keys))
  rb <- process_fit(fits_base[[key]], dat, B, seed = 5000 + match(key, keys))
  
  edagH_out[[length(edagH_out) + 1]] <- tibble(
    country = cc, sex = sx, year = rt$years,
    LE_dA = rt$edagH[, 1], edag_dA = rt$edagH[, 2], H_dA = rt$edagH[, 3])
  
  ## ---- improvement-series bands ----
  it  <- lapply(rt$draws, impr_mat, years = rt$years)
  itp <- lapply(rt$point, function(v) impr_mat(matrix(v, ncol = 1), rt$years))
  ib  <- lapply(rb$draws, impr_mat, years = rb$years)
  ibp <- lapply(rb$point, function(v) impr_mat(matrix(v, ncol = 1), rb$years))
  yy  <- it$dAPC$years
  
  for (s in c("dA", "dAC", "dAP", "dAPC")) {
    out[[length(out) + 1]] <- band_df(yy, itp[[s]]$impr[, 1], it[[s]]$impr, s, "Total", cc, sx)
  }
  ## baseline: only the full-model (dAPC) LE improvement is needed for the gap (H3/H5)
  out[[length(out) + 1]] <- band_df(yy, ibp$dAPC$impr[, 1], ib$dAPC$impr, "dAPC", "Baseline", cc, sx)
  out[[length(out) + 1]] <- band_df(yy,
                                    itp$dAPC$impr[, 1] - itp$dAP$impr[, 1], it$dAPC$impr - it$dAP$impr,
                                    "contrib_cohort", "Total", cc, sx)
  out[[length(out) + 1]] <- band_df(yy,
                                    itp$dAPC$impr[, 1] - itp$dAC$impr[, 1], it$dAPC$impr - it$dAC$impr,
                                    "contrib_period", "Total", cc, sx)
  out[[length(out) + 1]] <- band_df(yy,
                                    ibp$dAPC$impr[, 1] - itp$dAPC$impr[, 1], ib$dAPC$impr - it$dAPC$impr,
                                    "gap_dAPC", "Gap", cc, sx)
  
  ## ---- rolling-Shapley contributions on the improvement series (Total) ----
  ## Same two-player Shapley as the period table, but on the 10-yr rolling
  ## improvement. drift + period + cohort = dAPC improvement at every year.
  sh_drift_p  <- itp$dA$impr[, 1]
  sh_period_p <- 0.5 * ((itp$dAP$impr[, 1] - itp$dA$impr[, 1]) +
                          (itp$dAPC$impr[, 1] - itp$dAC$impr[, 1]))
  sh_cohort_p <- 0.5 * ((itp$dAC$impr[, 1] - itp$dA$impr[, 1]) +
                          (itp$dAPC$impr[, 1] - itp$dAP$impr[, 1]))
  sh_drift_d  <- it$dA$impr
  sh_period_d <- 0.5 * ((it$dAP$impr - it$dA$impr) + (it$dAPC$impr - it$dAC$impr))
  sh_cohort_d <- 0.5 * ((it$dAC$impr - it$dA$impr) + (it$dAPC$impr - it$dAP$impr))
  
  out[[length(out) + 1]] <- band_df(yy, sh_drift_p,  sh_drift_d,  "shap_drift",  "Total", cc, sx)
  out[[length(out) + 1]] <- band_df(yy, sh_period_p, sh_period_d, "shap_period", "Total", cc, sx)
  out[[length(out) + 1]] <- band_df(yy, sh_cohort_p, sh_cohort_d, "shap_cohort", "Total", cc, sx)
  
  ## ---- period change-decomposition (Total rates; reuses rt) ----
  for (pp in seq_len(nrow(periods))) {
    t0 <- periods$t0[pp]; t1 <- periods$t1[pp]; per <- periods$period[pp]
    
    dAd   <- period_change(rt$draws$dA,   rt$years, t0, t1)
    dAPd  <- period_change(rt$draws$dAP,  rt$years, t0, t1)
    dACd  <- period_change(rt$draws$dAC,  rt$years, t0, t1)
    dAPCd <- period_change(rt$draws$dAPC, rt$years, t0, t1)
    dAp   <- period_change(rt$point$dA,   rt$years, t0, t1)
    dAPp  <- period_change(rt$point$dAP,  rt$years, t0, t1)
    dACp  <- period_change(rt$point$dAC,  rt$years, t0, t1)
    dAPCp <- period_change(rt$point$dAPC, rt$years, t0, t1)
    
    ## Drift held fixed as the baseline (age+drift). Two-player Shapley over the
    ## non-linear period and cohort components: each = average of its marginal
    ## contribution across the two orderings. drift + period + cohort = total exactly.
    drift_d  <- dAd
    period_d <- 0.5 * ((dAPd - dAd) + (dAPCd - dACd))
    cohort_d <- 0.5 * ((dACd - dAd) + (dAPCd - dAPd))
    total_d  <- dAPCd
    drift_p  <- dAp
    period_p <- 0.5 * ((dAPp - dAp) + (dAPCp - dACp))
    cohort_p <- 0.5 * ((dACp - dAp) + (dAPCp - dAPp))
    total_p  <- dAPCp
    
    bw <- function(dv, pt, comp) {
      q <- quantile(dv, c(0.025, 0.975), na.rm = TRUE)
      tibble(country = cc, sex = sx, period = per, component = comp,
             estimate = pt, lower = q[1], upper = q[2])
    }
    
    decomp[[length(decomp) + 1]] <- bind_rows(
      bw(drift_d,  drift_p,  "drift"),
      bw(period_d, period_p, "period"),
      bw(cohort_d, cohort_p, "cohort"),
      bw(total_d,  total_p,  "total")
    )
    ## stash per-draw and per-point components for the 1983-2017 summary row
    decomp_draws[[length(decomp_draws) + 1]] <- tibble(
      country = cc, sex = sx, per = per,
      draw = seq_len(B),
      drift = drift_d, period = period_d, cohort = cohort_d, total = total_d
    )
    decomp_pts[[length(decomp_pts) + 1]] <- tibble(
      country = cc, sex = sx, per = per,
      drift = drift_p, period = period_p, cohort = cohort_p, total = total_p
    )
    
  }
}

## ---- full-window summary row (1983-2017) -----------------------------------
## Point = length-weighted MLE point (consistent with the sub-period rows);
## interval = MC quantiles of the length-weighted average computed WITHIN each
## draw, so the interval reflects the correlation across sub-periods.
dd  <- bind_rows(decomp_draws) %>% mutate(w = wlen[per])
pts <- bind_rows(decomp_pts)   %>% mutate(w = wlen[per])

summary_pts <- pts %>%
  group_by(country, sex) %>%
  summarise(across(c(drift, period, cohort, total), ~ sum(.x * w) / sum(w)),
            .groups = "drop") %>%
  tidyr::pivot_longer(c(drift, period, cohort, total),
                      names_to = "component", values_to = "estimate")

summary_long <- dd %>%
  group_by(country, sex, draw) %>%
  summarise(across(c(drift, period, cohort, total), ~ sum(.x * w) / sum(w)),
            .groups = "drop") %>%
  tidyr::pivot_longer(c(drift, period, cohort, total),
                      names_to = "component", values_to = "value") %>%
  group_by(country, sex, component) %>%
  summarise(lower = quantile(value, .025),
            upper = quantile(value, .975), .groups = "drop") %>%
  left_join(summary_pts, by = c("country", "sex", "component")) %>%
  mutate(period = "1983-2017") %>%
  dplyr::select(country, sex, period, component, estimate, lower, upper)

## additivity guard: drift + period + cohort = total on the summary point
chk <- summary_long %>%
  tidyr::pivot_wider(id_cols = c(country, sex),
                     names_from = component, values_from = estimate) %>%
  mutate(resid = drift + period + cohort - total)
stopifnot(max(abs(chk$resid)) < 1e-8)

## ---- write improvement-series bands ----------------------------------------
result <- bind_rows(out)
write.csv(result, file.path(models.dir, "LE_series_ci.csv"), row.names = FALSE)
cat("\nwrote", nrow(result), "rows to LE_series_ci.csv\n")

## ---- write change-decomposition (sub-periods + 1983-2017 summary) ----------
decomp_tbl <- bind_rows(bind_rows(decomp), summary_long) %>%
  mutate(component = factor(component,
                            levels = c("drift", "period", "cohort", "total")),
         period    = factor(period,
                            levels = c(periods$period, "1983-2017"))) %>%
  arrange(country, sex, period, component)
write.csv(decomp_tbl, file.path(models.dir, "LE_change_decomp.csv"), row.names = FALSE)

wide <- decomp_tbl %>%
  mutate(component = recode(component,
                            drift    = "Drift",
                            period   = "NL_period",
                            cohort   = "NL_cohort",
                            total    = "Total"),
         txt = sprintf("%.3f (%.3f, %.3f)", estimate, lower, upper)) %>%
  dplyr::select(country, sex, period, component, txt) %>%
  tidyr::pivot_wider(names_from = component, values_from = txt)

write.csv(wide, file.path(models.dir, "LE_change_decomp_wide.csv"), row.names = FALSE)

### save the e-dagger, Keyfitz-H results
edagH_tbl <- bind_rows(edagH_out)
write.csv(edagH_tbl, file.path(models.dir, "dA_edag_H.csv"), row.names = FALSE)

## supports the H1 sentence: H_dA should slope < 0 everywhere; edag_dA flat/slightly +
edagH_trend <- edagH_tbl %>% group_by(country, sex) %>%
  summarise(slope_edag = coef(lm(edag_dA ~ year))[2],
            slope_H    = coef(lm(H_dA    ~ year))[2], .groups = "drop")
cat(sprintf("H_dA slope<0 in %d/%d; edag_dA slope>=0 in %d/%d\n",
            sum(edagH_trend$slope_H < 0), nrow(edagH_trend),
            sum(edagH_trend$slope_edag >= 0), nrow(edagH_trend)))