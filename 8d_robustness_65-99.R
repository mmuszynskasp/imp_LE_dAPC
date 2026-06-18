###############################################################################
## ROBUSTNESS CHECK: partial LE 65-99 vs the main 65-94 specification.
###############################################################################
rm(list = ls())

library(HMDHFDplus)
library(Epi)
library(dplyr)
library(tidyr)
library(MASS)        
library(ggplot2)
## ---- settings --------------------------------------------------------------
out.dir   <- #"your/data/directory"  
plot.dir   <- #"your/plots/directory" 

username <- "your_HMD_username"        
password <- "your_HMD_password"          

## full study set, including the flagged GBRTENW
countries  <- c("CHE","DNK","FIN","FRATNP","GBRTENW","ITA","NLD","NOR","SWE")
years      <- 1980:2019
age_lo     <- 65
age_hi     <- 99            # single-year top age; the life table is closed here (99+)
B          <- 2000
TOL        <- 1e-2
ref_p      <- 1980 + 1/3

## cohort range widened for the higher top age and both Lexis triangles:
##   oldest: age 99 in 1980, lower triangle -> born 1880
##   youngest: age 65 in 2019, upper triangle -> born 1944
## (1880:1944 spans every triangle in the 65-99 x 1980-2019 box; the filter is
##  effectively a documented no-op but kept for clarity.)
ourcohorts <- 1880:1944

###############################################################################
## 1. Prepare Lexis data, single-year ages 65-99 -----------------------------
###############################################################################
## triangle centroids exactly as in 01_prepare_data.R
get_country <- function(cc) {
  exp <- readHMDweb(cc, "Exposures_lexis", username, password) %>%
    filter(Year %in% years, Age >= age_lo) %>%
    dplyr::select(-Total, -OpenInterval) %>%
    pivot_longer(c(Female, Male), names_to = "sex", values_to = "Exp")
  dth <- readHMDweb(cc, "Deaths_lexis", username, password) %>%
    filter(Year %in% years, Age >= age_lo) %>%
    dplyr::select(-Total, -OpenInterval) %>%
    pivot_longer(c(Female, Male), names_to = "sex", values_to = "Dx")
  
  raw <- exp %>%
    left_join(dth, by = c("Year", "Age", "Cohort", "sex")) %>%
    mutate(country = cc)
  
  ## single-year ages 65-99: keep triangles, restrict cohorts, shift centroids
  raw %>%
    filter(Age >= age_lo, Age <= age_hi, Cohort %in% ourcohorts) %>%
    mutate(older = ifelse(Cohort == Year - Age, 1L, 0L),
           Age   = ifelse(older == 1L, Age  + 1/3, Age  + 2/3),
           Year2 = ifelse(older == 1L, Year + 2/3, Year + 1/3)) %>%
    mutate(Year = Year2) %>% dplyr::select(-Year2) %>%
    mutate(Cohort = Year - Age)
}

lexis <- bind_rows(lapply(countries, get_country))

## defensive: no missing/zero exposure should enter the fit
if (any(is.na(lexis$Exp) | is.na(lexis$Dx)))
  warning("NA in Dx/Exp after assembly - check the left_join keys")
write.table(lexis, file.path(out.dir, "lexis_6599.csv"), sep = ",", row.names = FALSE)

###############################################################################
## 2. Fit the detrended APC model, ages 65-99 --------------------------------
###############################################################################
fits <- list()
for (cc in countries) for (sx in c("Female","Male")) {
  dat <- lexis %>%
    filter(country == cc, sex == sx) %>%
    transmute(A = Age, P = Year, D = round(Dx), Y = round(Exp)) %>%
    as.data.frame()
  m <- apc.fit(data = dat,
               npar = c(A = 10, P = 10, C = 10),                # same knot count as 02;
               model = "bs", parm = "APC",                      # apc.fit places knots over the
               dr.extr = "Y", ref.p = ref_p)                    # observed (now wider) age range
  fits[[paste(cc, sx, sep = "_")]] <- m
}
saveRDS(fits, file.path(out.dir, "fits_bs_6599.rds"))

###############################################################################
## 3. Decomposition helpers (identical to 05) ---------------------------------
###############################################################################
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
  Lx[K, ] <- lx[K, ] / M[K, ]          # open interval closed at the top age (99+)
  colSums(Lx)
}
make_drift_vec <- function(Xp, u, w) {
  Sw <- sum(w); Swu <- sum(w * u); Swuu <- sum(w * u * u)
  denom <- Sw * Swuu - Swu^2
  (Sw * colSums((w * u) * Xp) - Swu * colSums(w * Xp)) / denom
}
period_change <- function(x, years, t0, t1) {
  i0 <- match(t0, years); i1 <- match(t1, years)
  if (is.matrix(x)) (x[i1, ] - x[i0, ]) / (t1 - t0) else (x[i1] - x[i0]) / (t1 - t0)
}

periods <- data.frame(
  period = c("1983-1993","1993-2003","2003-2013","2013-2017"),
  t0 = c(1983,1993,2003,2013), t1 = c(1993,2003,2013,2017),
  stringsAsFactors = FALSE)

## process_fit: as in 05. All single-year ages 65-99 are model-derived; the
## life table is closed at age 99 (l99/m99) inside pLE_matrix. No 100+ group.
process_fit <- function(m, dat, B, seed) {
  b <- coef(m$Model); nm <- names(b)
  phi <- sum(residuals(m$Model, type = "pearson")^2) / m$Model$df.residual
  phi <- max(phi, 1)
  V   <- phi * vcov(m$Model)
  ia <- which(startsWith(nm, "MA")); ip <- which(startsWith(nm, "MPr")); ic <- which(startsWith(nm, "MCr"))
  X  <- model.matrix(m$Model)
  Xa <- X[, ia, drop=FALSE]; Xp <- X[, ip, drop=FALSE]; Xc <- X[, ic, drop=FALSE]
  refP <- as.numeric(m$Ref["Per"]); u <- dat$P - refP; w <- dat$Y
  dvec <- make_drift_vec(Xp, u, w)
  drift_computed <- exp(sum(dvec * b[ip]))
  if (abs(drift_computed - as.numeric(m$Drift[1,1])) / as.numeric(m$Drift[1,1]) > TOL)
    stop("drift mismatch")
  
  ageint <- floor(dat$A); yearint <- floor(dat$P)
  yrs <- sort(unique(yearint)); ags <- sort(unique(ageint))   # 65..99 (35 ages)
  cellf <- factor(paste(yearint, ageint, sep = "_"))
  denom <- rowsum(w, cellf)[, 1]
  if (any(denom <= 0)) stop("zero exposure in a (year,age) cell - cannot form a rate")
  
  le_from_lr <- function(LR) {
    rate  <- exp(LR)
    numer <- rowsum(rate * w, cellf)
    cellrate <- numer / denom
    LE <- matrix(NA_real_, length(yrs), ncol(LR))
    for (iy in seq_along(yrs)) {
      rn <- paste(yrs[iy], ags, sep = "_")
      mx_vec <- cellrate[rn, , drop = FALSE]              # 65..99 rates (35 ages)
      if (anyNA(mx_vec)) stop("missing (year,age) cell in cellrate")
      LE[iy, ] <- pLE_matrix(mx_vec)                      # closed at 99+
    }
    LE
  }
  
  set.seed(seed)
  Bd <- MASS::mvrnorm(B, b, V)
  Ea <- Xa %*% t(Bd[, ia, drop=FALSE]); Ep <- Xp %*% t(Bd[, ip, drop=FALSE]); Ec <- Xc %*% t(Bd[, ic, drop=FALSE])
  ld <- as.vector(Bd[, ip, drop=FALSE] %*% dvec); U <- outer(u, ld)
  draws <- list(dAPC = le_from_lr(Ea+Ep+Ec), dAP = le_from_lr(Ea+Ep),
                dAC  = le_from_lr(Ea+U+Ec),  dA  = le_from_lr(Ea+U))
  ea <- Xa %*% b[ia]; ep <- Xp %*% b[ip]; ec <- Xc %*% b[ic]; um <- u * sum(dvec * b[ip])
  point <- list(dAPC = le_from_lr(matrix(ea+ep+ec,ncol=1))[,1], dAP = le_from_lr(matrix(ea+ep,ncol=1))[,1],
                dAC  = le_from_lr(matrix(ea+um+ec,ncol=1))[,1], dA = le_from_lr(matrix(ea+um,ncol=1))[,1])
  list(years = yrs, draws = draws, point = point)
}

###############################################################################
## 4. Decompose each country x sex, build the comparison table ----------------
###############################################################################
decomp <- list()
for (key in names(fits)) {
  parts <- strsplit(key, "_")[[1]]; cc <- parts[1]; sx <- parts[2]
  
  dat <- lexis %>% filter(country == cc, sex == sx) %>%
    transmute(A = Age, P = Year, D = round(Dx), Y = round(Exp)) %>% as.data.frame()
  
  rt <- process_fit(fits[[key]], dat, B, seed = 7000 + match(key, names(fits)))
  
  for (pp in seq_len(nrow(periods))) {
    t0 <- periods$t0[pp]; t1 <- periods$t1[pp]; per <- periods$period[pp]
    dAd<-period_change(rt$draws$dA,rt$years,t0,t1);  dAPd<-period_change(rt$draws$dAP,rt$years,t0,t1)
    dACd<-period_change(rt$draws$dAC,rt$years,t0,t1);dAPCd<-period_change(rt$draws$dAPC,rt$years,t0,t1)
    dAp<-period_change(rt$point$dA,rt$years,t0,t1);  dAPp<-period_change(rt$point$dAP,rt$years,t0,t1)
    dACp<-period_change(rt$point$dAC,rt$years,t0,t1);dAPCp<-period_change(rt$point$dAPC,rt$years,t0,t1)
    period_d <- 0.5*((dAPd-dAd)+(dAPCd-dACd)); cohort_d <- 0.5*((dACd-dAd)+(dAPCd-dAPd))
    period_p <- 0.5*((dAPp-dAp)+(dAPCp-dACp)); cohort_p <- 0.5*((dACp-dAp)+(dAPCp-dAPp))
    bw <- function(dv, pt, comp) {
      q <- quantile(dv, c(0.025,0.975), na.rm = TRUE)
      tibble(country=cc, sex=sx, period=per, component=comp, estimate=pt, lower=q[1], upper=q[2])
    }
    decomp[[length(decomp)+1]] <- bind_rows(
      bw(dAd, dAp, "drift"), bw(period_d, period_p, "period"),
      bw(cohort_d, cohort_p, "cohort"), bw(dAPCd, dAPCp, "total"))
  }
}

decomp_6599 <- bind_rows(decomp) %>%
  mutate(age_range = "65-99",
         component = factor(component, levels = c("drift","period","cohort","total")),
         period    = factor(period, levels = periods$period)) %>%
  arrange(country, sex, period, component)

write.csv(decomp_6599, file.path(out.dir, "LE_change_decomp_6599.csv"), row.names = FALSE)

###############################################################################
## 5. Side-by-side with the main 65-94 decomposition --------------------------
##    (reads LE_change_decomp.csv from 05; joins the point estimates + CIs)
###############################################################################
main94 <- read.csv(file.path(out.dir, "LE_change_decomp.csv")) %>%
  mutate(age_range = "65-94") %>%
  dplyr::select(country, sex, period, component,
                est_94 = estimate, lo_94 = lower, hi_94 = upper)

compare <- decomp_6599 %>%
  dplyr::select(country, sex, period, component,
                est_99 = estimate, lo_99 = lower, hi_99 = upper) %>%
  mutate(period = as.character(period), component = as.character(component)) %>%
  left_join(main94 %>% mutate(period = as.character(period),
                              component = as.character(component)),
            by = c("country","sex","period","component")) %>%
  mutate(diff = est_99 - est_94,
         ## do the 94% intervals overlap?
         ci_overlap = !(lo_99 > hi_94 | lo_94 > hi_99)) %>%
  arrange(country, sex, period, component)

## guard: every 65-99 row must have matched a 65-94 row (label mismatches -> NA)
if (anyNA(compare$est_94))
  stop("unmatched rows in the 65-94 vs 65-99 join - check period/component labels")

write.csv(compare, file.path(out.dir, "decomp_94_vs_6599.csv"), row.names = FALSE)

###############################################################################
## 6. Quick read: is the decomposition stable to the cap? ---------------------
###############################################################################
cat("\nmax |65-99  minus  65-94| difference, by component:\n")
print(compare %>% group_by(component) %>%
        summarise(max_abs_diff = max(abs(diff)),
                  median_abs_diff = median(abs(diff)),
                  ci_overlap_rate = mean(ci_overlap), .groups = "drop") %>%
        as.data.frame(), row.names = FALSE)

cat("\nsign agreement of point estimates (period & cohort - qualitative story):\n")
print(compare %>% filter(component %in% c("period","cohort")) %>%
        summarise(n = n(),
                  same_sign = mean(sign(est_99) == sign(est_94)),
                  .groups = "drop") %>% as.data.frame(), row.names = FALSE)

cat("\nlargest 12 absolute differences:\n")
print(compare %>% arrange(desc(abs(diff))) %>%
        dplyr::select(country, sex, period, component, est_94, est_99, diff, ci_overlap) %>%
        head(12) %>% as.data.frame(), row.names = FALSE)

###############################################################################
## Robustness figure: per-cell decomposition under the main 65-94 spec vs the
## refit 65-99 spec. Each point is one country x sex x period x component cell;
## the dashed line is the identity (y = x). Points on the line are unchanged by
## the cap; systematic departure (the drift panel) is the modest attenuation.
##
## Reads decomp_94_vs_6599.csv written by robustness_6599.R.
###############################################################################

## component order + readable facet labels (period before cohort, as in the paper)
comp_levels <- c("drift", "period", "cohort", "total")
comp_labels <- c(drift  = "Drift",
                 period = "Non-linear period",
                 cohort = "Non-linear cohort",
                 total  = "Total (dAPC)")

dat <- compare %>%
  ## drop "total" here if you prefer a 3-panel figure of the additive parts:
  # filter(component != "total") %>%
  mutate(component = factor(component, levels = comp_levels),
         sex       = factor(sex, levels = c("Female", "Male")))

## optional: label the largest departures from the identity line (needs ggrepel)
## colourblind-safe (Okabe-Ito) for the two sexes
sex_cols <- c(Female = "#E69F00", Male = "#0072B2")

p <- ggplot(dat, aes(est_94, est_99)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey55", linewidth = 0.4) +
  geom_point(aes(colour = sex, shape = sex), size = 1.9, alpha = 0.85) +
  facet_wrap(~ component, scales = "free",
             labeller = as_labeller(comp_labels)) +
  scale_colour_manual(values = sex_cols, name = NULL) +
  scale_shape_manual(values = c(Female = 16, Male = 17), name = NULL) +
  labs(x = "Average annual change, ages 65-94 (years per year)",
       y = "Average annual change, ages 65-99 (years per year)") +
  theme_bw(base_size = 10) +
  theme(aspect.ratio   = 1,                 # square panels; identity line is exact
        panel.grid.minor = element_blank(),
        legend.position  = "top",
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text       = element_text(face = "bold"))


ggsave(file.path(plot.dir, "fig_robustness_6599_scatter.pdf"),
       p, width = 7.0, height = 7.0)
ggsave(file.path(plot.dir, "fig_robustness_6599_scatter.png"),
       p, width = 7.0, height = 7.0, dpi = 320)
