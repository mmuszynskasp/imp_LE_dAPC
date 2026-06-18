###############################################################################
## Project:  Why Life Expectancy Improvements Slowed (dAPC / LE)
## Aim:      Monte Carlo confidence bands for the LE-improvement series AND the
##           period change-decomposition, from one consistent run.
##           TOTAL RATES ONLY (baseline / "indiv" estimation removed).
##           Series, per country x sex (LE_series_ci_1960.csv):
##             dA, dAC, dAP, dAPC          (Total rates)
##             contrib_cohort = dAPC - dAP (Total; non-linear cohort)
##             contrib_period = dAPC - dAC (Total; non-linear period)
##             shap_drift, shap_period, shap_cohort (Total; rolling Shapley, sum to dAPC)
##           Change-decomposition over fixed periods (LE_change_decomp_1960.csv):
##             average annual change (LE(t1)-LE(t0))/(t1-t0), Total rates,
##             split into drift / period / cohort / total, with CIs.
## Inputs:   <models.dir>/fits_bs1960.rds   (total, from 02)
##           <data.dir>/lexis_1960.csv      (to map model rows to A, P, Y)
## Outputs:  <models.dir>/LE_series_ci_1960.csv
##           <models.dir>/LE_change_decomp_1960.csv  (+ _wide_1960.csv)
## Depends:  dplyr, tidyr, MASS, Epi (for the fits)
###############################################################################
rm(list = ls())

library(dplyr)
library(tidyr)
library(MASS)
library(Epi)

data.dir   <- #"your/data/directory"  
models.dir <- #"your/output/directory"  


B    <- 2000          # number of draws
TOL  <- 1e-2          # relative tolerance for the drift validation

fits_tot  <- readRDS(file.path(models.dir, "fits_bs1960.rds"))
lexis     <- read.table(file.path(data.dir, "lexis_1960.csv"), sep = ",", header = TRUE)

## periods for the change-decomposition (Total rates). The last is shorter (4y)
## -> its CI is wider; that is expected, not instability.
periods <- data.frame(
  period = c("1963-1973","1973-1983", "1983-1993", "1993-2003", "2003-2013", "2013-2017"),
  t0     = c(1963,1973, 1983, 1993, 2003, 2013),
  t1     = c(1973,1983,1993, 2003, 2013, 2017),
  stringsAsFactors = FALSE
)

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
  list(years = yrs, draws = draws, point = point,
       drift = c(computed = drift_computed, reported = drift_reported))
}

## ---- main loop over country x sex ------------------------------------------
keys   <- names(fits_tot)
out    <- list()    # improvement-series bands
decomp <- list()    # period change-decomposition

for (key in keys) {
  parts <- strsplit(key, "_")[[1]]; cc <- parts[1]; sx <- parts[2]
  dat <- lexis %>% filter(country == cc, sex == sx) %>%
    transmute(A = Age, P = Year, D = round(Dx), Y = round(Exp)) %>% as.data.frame()
  
  rt <- process_fit(fits_tot[[key]], dat, B, seed = 1000 + match(key, keys))
  cat(sprintf("%-16s drift Total %.4f (rep %.4f)\n",
              key, rt$drift[1], rt$drift[2]))
  
  ## ---- improvement-series bands (Total only) ----
  it  <- lapply(rt$draws, impr_mat, years = rt$years)
  itp <- lapply(rt$point, function(v) impr_mat(matrix(v, ncol = 1), rt$years))
  yy  <- it$dAPC$years
  
  for (s in c("dA", "dAC", "dAP", "dAPC")) {
    out[[length(out) + 1]] <- band_df(yy, itp[[s]]$impr[, 1], it[[s]]$impr, s, "Total", cc, sx)
  }
  out[[length(out) + 1]] <- band_df(yy,
                                    itp$dAPC$impr[, 1] - itp$dAP$impr[, 1], it$dAPC$impr - it$dAP$impr,
                                    "contrib_cohort", "Total", cc, sx)
  out[[length(out) + 1]] <- band_df(yy,
                                    itp$dAPC$impr[, 1] - itp$dAC$impr[, 1], it$dAPC$impr - it$dAC$impr,
                                    "contrib_period", "Total", cc, sx)
  
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
  }
}

## ---- write improvement-series bands ----------------------------------------
result <- bind_rows(out)
write.csv(result, file.path(models.dir, "LE_series_ci_1960.csv"), row.names = FALSE)

## ---- write change-decomposition (tidy + manuscript-ready wide) -------------
decomp_tbl <- bind_rows(decomp) %>%
  mutate(component = factor(component,
                            levels = c("drift", "period", "cohort", "total")),
         period    = factor(period, levels = periods$period)) %>%
  arrange(country, sex, period, component)
write.csv(decomp_tbl, file.path(models.dir, "LE_change_decomp_1960.csv"), row.names = FALSE)

wide <- decomp_tbl %>%
  mutate(component = recode(component,
                            drift    = "Drift",
                            period   = "NL_period",
                            cohort   = "NL_cohort",
                            total    = "Total"),
         txt = sprintf("%.3f (%.3f, %.3f)", estimate, lower, upper)) %>%
  dplyr::select(country, sex, period, component, txt) %>%
  tidyr::pivot_wider(names_from = component, values_from = txt)

write.csv(wide, file.path(models.dir, "LE_change_decomp_wide_1960.csv"), row.names = FALSE)