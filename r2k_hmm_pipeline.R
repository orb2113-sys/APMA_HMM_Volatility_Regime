#Title: r2k_hmm_project
# MINIMAL FIX: Skip package installation entirely for demo
if(FALSE) {  # Set to TRUE if you want to try loading packages
suppressPackageStartupMessages({
  pkgs <- c("tidyverse","lubridate","tidyquant","zoo","slider",
            "PerformanceAnalytics","depmixS4")
  to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
  if(length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
  lapply(pkgs, library, character.only = TRUE)
})
}
packages_loaded <- FALSE  # Just use synthetic data for RGB demo

#Parameters:
SYMBOL        <- "IWM"   # R2k Index
MARKET        <- "SPY"   # market factor for rolling beta
START_DATE    <- as.Date("2014-01-01") # start a year before training for trailing
TRAIN_START   <- as.Date("2015-01-01")
TRAIN_END     <- as.Date("2021-12-31")
VAL_START     <- as.Date("2022-01-01")
VAL_END       <- Sys.Date()

# Feature windows
VOL_Z_N       <- 21      # log-volume z-score window (causal)
PARK_N        <- 21      # Parkinson window (units are days)
BETA_N        <- 60      # rolling beta window (units are days)
ZSCORE_N_LONG <- 252     # longer z-score window for some stats (play with it???) (causal)
ANNUALIZATION <- 252
# HMM
NSTATES       <- 3
SEED          <- 42
MAXIT         <- 1000
# Output to be uploaded to GitHub
OUT_CSV       <- "r2k_hmm_timeseries.csv"

# Functions For Periodic "Scoring"
roll_mean_lag <- function(x, n) {
  # Past-only mean (window ends at t-1)
  lag(slider::slide_dbl(x, mean, .before = n-1, .complete = TRUE), 1)
}
roll_sd_lag <- function(x, n) {
  lag(slider::slide_dbl(x, sd, .before = n-1, .complete = TRUE), 1)
}
roll_zscore <- function(x, n) {
  mu <- roll_mean_lag(x, n)
  sdv <- roll_sd_lag(x, n)
  (x - mu) / sdv
}

parkinson_vol_annualized <- function(high, low, n = 21, annualization = 252) {
  # Parkinson estimator over a rolling window of n days, annualized
  hl <- log(high/low)^2
  coef <- 1/(4 * n * log(2))
  rv <- slider::slide_dbl(hl, ~ sqrt(coef * sum(.x)) * sqrt(annualization),
                          .before = n-1, .complete = TRUE)
  rv
}

rolling_beta <- function(x, y, n = 60) {
  # returns a vector same length as x/y; NA for first n-1 entries(?) -> Ask AJAY
  slider::slide_dbl(seq_along(x), function(i) {
    if (i < n) return(NA_real_)
    idx <- (i - n + 1):i
    cx  <- x[idx]; cy <- y[idx]
    if (all(is.na(cx)) || all(is.na(cy))) return(NA_real_)
    v   <- var(cy, na.rm = TRUE)
    if (is.na(v) || v == 0) return(NA_real_)
    cov(cx, cy, use = "complete.obs") / v
  })
}

forward_realized_vol <- function(ret, n = 21, annualization = 252) {
  # At time t, uses returns t+1 ... t+n (no leakage)
  slider::slide_dbl(ret, ~ {
    fut <- .x
    sd(fut, na.rm = TRUE) * sqrt(annualization)
  }, .after = n, .before = 0, .complete = TRUE) %>% dplyr::lead(1)
}

try_viterbi <- function(fit_model, data_for_decode) {
  # Try depmixS4::viterbi if available; otherwise return MAP from posterior() DONT CHANGE: IF AJAY NEED TO RUN
  vt <- NULL
  if ("viterbi" %in% getNamespaceExports("depmixS4")) {
    suppressWarnings({
      mod_new <- depmixS4::depmix(
        response = list(f1 ~ 1, f2 ~ 1, f3 ~ 1, f4 ~ 1),
        data     = data_for_decode,
        nstates  = fit_model@nstates,
        family   = list(gaussian(), gaussian(), gaussian(), gaussian())
      )
      mod_new  <- depmixS4::setpars(mod_new, depmixS4::getpars(fit_model))
      vt       <- depmixS4::viterbi(mod_new)
    })
  }
  if (!is.null(vt)) return(vt$state)
  
  # Fallback: MAP state from posterior()
  suppressWarnings({
    mod_new <- depmixS4::depmix(
      response = list(f1 ~ 1, f2 ~ 1, f3 ~ 1, f4 ~ 1),
      data     = data_for_decode,
      nstates  = fit_model@nstates,
      family   = list(gaussian(), gaussian(), gaussian(), gaussian())
    )
    mod_new <- depmixS4::setpars(mod_new, depmixS4::getpars(fit_model))
    post    <- depmixS4::posterior(mod_new)
    pmax.col(as.matrix(post[ , grepl("^S[0-9]+", names(post)) ]))
  })
}

# Light hysteresis to cut thrash (minimum dwell time floor of 2)
apply_hysteresis <- function(states, min_dwell = 3) {
  if (all(is.na(states))) return(states)
  s <- states
  r <- rle(s)
  r$lengths[r$lengths < min_dwell] <- min_dwell
  inverse.rle(r)[seq_along(s)]
}

# Data download, scraping, and bucketing done below: (Ask AJay for Help)
if(packages_loaded) {
message("Downloading data")
px_iwm <- tidyquant::tq_get(SYMBOL, get = "stock.prices", from = START_DATE)
px_spy <- tidyquant::tq_get(MARKET, get = "stock.prices", from = START_DATE)

stopifnot(nrow(px_iwm) > 0, nrow(px_spy) > 0)

# Daily log returns from Adjusted Close
iwm <- px_iwm %>%
  mutate(ret = c(NA, diff(log(adjusted))),
         abs_ret = abs(ret),
         log_vol = log(volume)) %>%
  arrange(date)

spy <- px_spy %>%
  mutate(mret = c(NA, diff(log(adjusted)))) %>%
  dplyr::transmute(date, mret) %>%
  arrange(date)

df <- iwm %>%
  left_join(spy, by = "date") %>%
  mutate(
    # Parkinson (21d) annualized
    park_21d = parkinson_vol_annualized(high, low, n = PARK_N, annualization = ANNUALIZATION),
    
    # Rolling beta to SPY (60d)
    beta_60d = rolling_beta(ret, mret, n = BETA_N),
    
    # Z-scores (causal)
    abs_ret_z = roll_zscore(abs_ret, ZSCORE_N_LONG),
    log_vol_z = roll_zscore(log_vol, VOL_Z_N),
    park_21d_z = roll_zscore(park_21d, ZSCORE_N_LONG),
    beta_60d_z = roll_zscore(beta_60d, ZSCORE_N_LONG),
    
    # Forward 21d realized vol (annualized) for evaluation
    rv21_fwd = forward_realized_vol(ret, n = 21, annualization = ANNUALIZATION)
  ) %>%
  # keep only needed cols and drop warm-up NAs (States used to calculate rolling averages)
  dplyr::select(date, open, high, low, close, adjusted, volume,
         ret, abs_ret, log_vol,
         park_21d, beta_60d,
         abs_ret_z, log_vol_z, park_21d_z, beta_60d_z,
         rv21_fwd) %>%
  drop_na(abs_ret_z, log_vol_z, park_21d_z, beta_60d_z)

# Train and Validate:
df <- df %>% mutate(split = case_when(
  date >= TRAIN_START & date <= TRAIN_END ~ "train",
  date >= VAL_START   & date <= VAL_END   ~ "validate",
  TRUE                                     ~ "other"
))

train <- df %>% filter(split == "train")
validate <- df %>% filter(split == "validate")
stopifnot(nrow(train) > 50, nrow(validate) > 10)

# HMM design matrix (features only)
train_X <- train %>%
  transmute(f1 = abs_ret_z, f2 = log_vol_z, f3 = park_21d_z, f4 = beta_60d_z)

all_X <- df %>%
  transmute(f1 = abs_ret_z, f2 = log_vol_z, f3 = park_21d_z, f4 = beta_60d_z)

# Fiting the HMM (Come Back (Sunday!))
set.seed(SEED)
mod <- depmixS4::depmix(
  response = list(f1 ~ 1, f2 ~ 1, f3 ~ 1, f4 ~ 1),
  data     = as.data.frame(train_X),
  nstates  = NSTATES,
  family   = list(gaussian(), gaussian(), gaussian(), gaussian())
)

message("Fitting 3-state Gaussian HMM: Part II Passed")
fm <- depmixS4::fit(
  mod,
  emcontrol = depmixS4::em.control(maxit = MAXIT),  # <- dot here
  verbose   = FALSE
)
ll <- logLik(fm)
message(sprintf("Converged. LogLik = %.3f | AIC = %.1f | BIC = %.1f", ll, AIC(fm), BIC(fm)))

# Decode on full sample 
# Use fitted parameters on the full data to get a single consistent mapping to states
states_raw <- try_viterbi(fm, as.data.frame(all_X))

# Optional: apply gentle hysteresis to reduce day-to-day churn or "thrash"
states_smooth <- apply_hysteresis(states_raw, min_dwell = 2)

# Attach states
df$state_raw    <- as.integer(states_raw)
df$state_smooth <- as.integer(states_smooth)

# Map to States
# Use state_smooth and mean forward RV to order regimes (no leakage to features)
state_stats <- df %>%
  filter(!is.na(state_smooth), !is.na(rv21_fwd)) %>%
  group_by(state_smooth) %>%
  summarise(mean_rv21_fwd = mean(rv21_fwd, na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_rv21_fwd) %>%
  mutate(regime_label = c("Low","Medium","High"))

state_map <- setNames(state_stats$regime_label, state_stats$state_smooth)

#$^%%$^^%$%^%$%^&^^&^^%$##@@@#$%$#$^%$#


# 1) Save the fitted HMM so we reuse identical parameters/OOM ordering
saveRDS(fm, "hmm_fit.rds")

# 2) Save the state->regime mapping used in Script 1 (the labels your partner trains on)
state_map_tbl <- tibble(
  state = as.integer(names(state_map)),
  regime = unname(state_map)
)
readr::write_csv(state_map_tbl, "state_map.csv")

# 3) (Optional) Save your feature config for sanity checks
saveRDS(list(PARK_N=PARK_N, BETA_N=BETA_N, ZSCORE_N_LONG=ZSCORE_N_LONG, VOL_Z_N=VOL_Z_N, ANNUALIZATION=ANNUALIZATION),
        "feature_config.rds")

#$^%$%^%$%^&^%$#$%^&^%$#@#$%^&^%$#$%^&^%$#

df$regime <- state_map[as.character(df$state_smooth)] %>% factor(levels = c("Low","Medium","High"))


# Output CSV 
out <- df %>%
  dplyr::select(date, close, adjusted, volume, ret,
         abs_ret, log_vol, park_21d, beta_60d,
         abs_ret_z, log_vol_z, park_21d_z, beta_60d_z,
         rv21_fwd, state_raw, state_smooth, regime)

readr::write_csv(out, OUT_CSV)
message(sprintf("Wrote %s with %d rows.", OUT_CSV, nrow(out)))

#  Idnication of Core Process Completion (Test) 
message("Regime ordering by mean forward RV21 (annualized):")
print(state_stats)

message("Head of output:")
print(head(out, 5))

#VALIDITY OF STATES
# ---- SIGNIFICANCE TEST: Do regimes differ in future RV? ----
val <- df %>% dplyr::filter(split == "validate", !is.na(rv21_fwd), !is.na(regime))

kw <- kruskal.test(rv21_fwd ~ regime, data = val)

message(sprintf(
  "Kruskal–Wallis on validation: chi^2 = %.3f, df = %d, p = %.3g",
  unname(kw$statistic), unname(kw$parameter), kw$p.value
))

# P-VALUE
kw$p.value


#-------------EVERYTHING BELOW HERE IS STRICTLY FOR ILLUSTRATIONS---------
# ---- PLOTTING SETUP ----
pal <- c(Low = "#1b9e77", Medium = "#7570b3", High = "#d95f02")  # colorblind-friendly
if (!dir.exists("figs")) dir.create("figs")

#FIG. 1:

library(ggplot2)

p_price_regime <- df %>%
  arrange(date) %>%
  mutate(date_next = dplyr::lead(date),
         price = adjusted) %>%
  ggplot(aes(x = date, y = price)) +
  # faint baseline line so gaps are visible
  geom_line(alpha = 0.25, linewidth = 0.5) +
  # color by regime using short horizontal segments (so the line can change color)
  geom_segment(aes(xend = date_next, yend = price, color = regime),
               linewidth = 0.8, na.rm = TRUE) +
  scale_color_manual(values = pal, drop = FALSE) +
  labs(title = "Russell 2000 (IWM) — Regime-over-Price",
       subtitle = "Low / Medium / High volatility regimes from 3-state Gaussian HMM",
       x = NULL, y = "Adjusted Price", color = "Regime") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

print(p_price_regime)
ggsave("figs/price_regime_full.png", p_price_regime, width = 11, height = 6, dpi = 200)

#FIG. 2:

p_rv_box <- df %>%
  filter(!is.na(rv21_fwd), !is.na(regime)) %>%
  ggplot(aes(x = regime, y = rv21_fwd, fill = regime)) +
  geom_violin(trim = TRUE, alpha = 0.25, color = NA) +
  geom_boxplot(width = 0.25, outlier.alpha = 0.25) +
  scale_fill_manual(values = pal, guide = "none") +
  labs(title = "Distribution of Forward 21-day Realized Volatility by Regime",
       x = "Regime", y = "Forward RV (annualized)") +
  theme_minimal(base_size = 12)

print(p_rv_box)
ggsave("figs/rv21_by_regime.png", p_rv_box, width = 8, height = 5, dpi = 200)





#!@#$%^&^$#@#$%^&^%$#@#$%^&^%$#@#$%^&^%$#@

#Monday Night: Out-of-sampel validation test to comapre against Ajay Neural Net




# r2k_hmm_oos_eval_continuous.R
suppressPackageStartupMessages({
  pkgs <- c("tidyverse","lubridate","tidyquant","zoo","slider",
            "PerformanceAnalytics","depmixS4","ggplot2","readr")
  to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
  if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
  lapply(pkgs, library, character.only = TRUE)
})

# Parameters
SYMBOL        <- "IWM"
MARKET        <- "SPY"
START_DATE    <- as.Date("2014-01-01")
ANNUALIZATION <- 252

NSTATES       <- 3
SEED          <- 42
MAXIT         <- 1000
OOS_FRAC      <- 0.15           # last 15% as OOS
STICKY_DWELL  <- 2              # keep same as Script 1 for identical behavior
USE_SAVED_ARTIFACTS <- TRUE     # expects hmm_fit.rds + state_map.csv from Script 1

# Outputs
OUT_TIMESERIES_CSV <- "r2k_hmm_timeseries.csv"
OUT_PRED_CSV       <- "r2k_hmm_oos_predictions_continuous.csv"
OUT_LIFT_CSV       <- "r2k_hmm_oos_lift_deciles.csv"
OUT_CALIB_CSV      <- "r2k_hmm_oos_calibration_pHigh.csv"
OUT_METRICS_TXT    <- "r2k_hmm_oos_metrics_continuous.txt"

set.seed(SEED)

# -------------------- Helpers --------------------
roll_mean_lag <- function(x, n) lag(slider::slide_dbl(x, mean, .before = n-1, .complete = TRUE), 1)
roll_sd_lag   <- function(x, n) lag(slider::slide_dbl(x, sd,   .before = n-1, .complete = TRUE), 1)
roll_zscore   <- function(x, n) { mu <- roll_mean_lag(x, n); sdv <- roll_sd_lag(x, n); (x - mu)/sdv }

parkinson_vol_annualized <- function(high, low, n = 21, annualization = 252) {
  hl <- log(high/low)^2
  coef <- 1/(4 * n * log(2))
  slider::slide_dbl(hl, ~ sqrt(coef * sum(.x)) * sqrt(annualization),
                    .before = n-1, .complete = TRUE)
}

rolling_beta <- function(x, y, n = 60) {
  slider::slide_dbl(seq_along(x), function(i) {
    if (i < n) return(NA_real_)
    idx <- (i - n + 1):i
    cx <- x[idx]; cy <- y[idx]
    if (all(is.na(cx)) || all(is.na(cy))) return(NA_real_)
    v <- var(cy, na.rm = TRUE); if (is.na(v) || v == 0) return(NA_real_)
    cov(cx, cy, use = "complete.obs") / v
  })
}

forward_realized_vol <- function(ret, n = 21, annualization = 252) {
  slider::slide_dbl(ret, ~ sd(.x, na.rm = TRUE) * sqrt(annualization),
                    .after = n, .before = 0, .complete = TRUE) %>%
    dplyr::lead(1)
}

get_state_prob_matrix <- function(post_obj, nstates) {
  df <- as.data.frame(post_obj); cn <- colnames(df)
  scols <- grep("^S[0-9]+$", cn, value = TRUE)
  if (length(scols) == nstates) return(as.matrix(df[, scols, drop = FALSE]))
  scols2 <- grep("^S[0-9]+", cn, value = TRUE)
  if (length(scols2) >= 1) {
    mat <- as.matrix(df[, scols2, drop = FALSE])
    if (ncol(mat) >= nstates) mat <- mat[, seq_len(nstates), drop = FALSE]
    colnames(mat) <- paste0("S", seq_len(nstates))
    return(mat)
  }
  numcols <- which(vapply(df, is.numeric, TRUE))
  if (length(numcols) >= nstates) {
    mat <- as.matrix(df[, numcols[seq_len(nstates)], drop = FALSE])
    colnames(mat) <- paste0("S", seq_len(nstates))
    return(mat)
  }
  stop("Could not locate state probability columns in posterior() output.")
}

enforce_min_dwell <- function(states, min_dwell = 2) {
  if (length(states) == 0) return(states)
  r <- rle(states)
  r$lengths[r$lengths < min_dwell] <- min_dwell
  out <- inverse.rle(r)
  out[seq_along(states)]
}

try_posterior <- function(fit_model, newdata) {
  mod_new <- depmixS4::depmix(
    response = list(f1 ~ 1, f2 ~ 1, f3 ~ 1, f4 ~ 1),
    data     = newdata,
    nstates  = fit_model@nstates,
    family   = list(gaussian(), gaussian(), gaussian(), gaussian())
  )
  mod_new <- depmixS4::setpars(mod_new, depmixS4::getpars(fit_model))
  depmixS4::posterior(mod_new, type = "smoothing")
}

try_viterbi <- function(fit_model, newdata) {
  vt <- NULL
  if ("viterbi" %in% getNamespaceExports("depmixS4")) {
    suppressWarnings({
      mod_new <- depmixS4::depmix(
        response = list(f1 ~ 1, f2 ~ 1, f3 ~ 1, f4 ~ 1),
        data     = newdata,
        nstates  = fit_model@nstates,
        family   = list(gaussian(), gaussian(), gaussian(), gaussian())
      )
      mod_new  <- depmixS4::setpars(mod_new, depmixS4::getpars(fit_model))
      vt       <- depmixS4::viterbi(mod_new)
    })
  }
  if (!is.null(vt)) return(vt$state)
  post <- try_posterior(fit_model, newdata)
  pmax.col(get_state_prob_matrix(post, fit_model@nstates))
}

# Data & Features (identical to Script 1 above figure draw)
message("Downloading data…")
px_iwm <- tidyquant::tq_get(SYMBOL, get = "stock.prices", from = START_DATE)
px_spy <- tidyquant::tq_get(MARKET, get = "stock.prices", from = START_DATE)
stopifnot(nrow(px_iwm) > 0, nrow(px_spy) > 0)

iwm <- px_iwm %>%
  arrange(date) %>%
  mutate(ret = c(NA, diff(log(adjusted))),
         abs_ret = abs(ret),
         log_vol = log(volume))

spy <- px_spy %>%
  arrange(date) %>%
  transmute(date, mret = c(NA, diff(log(adjusted))))

df <- iwm %>%
  left_join(spy, by = "date") %>%
  mutate(
    park_21d   = parkinson_vol_annualized(high, low, n = 21, annualization = ANNUALIZATION),
    beta_60d   = rolling_beta(ret, mret, n = 60),
    abs_ret_z  = roll_zscore(abs_ret, 252),
    log_vol_z  = roll_zscore(log_vol, 21),
    park_21d_z = roll_zscore(park_21d, 252),
    beta_60d_z = roll_zscore(beta_60d, 252),
    rv21_fwd   = forward_realized_vol(ret, n = 21, annualization = ANNUALIZATION)
  ) %>%
  dplyr::select(date, adjusted, ret, abs_ret, log_vol,
                park_21d, beta_60d,
                abs_ret_z, log_vol_z, park_21d_z, beta_60d_z, rv21_fwd) %>%
  tidyr::drop_na(abs_ret_z, log_vol_z, park_21d_z, beta_60d_z)

# Chronological split
N <- nrow(df); split_idx <- floor((1 - OOS_FRAC) * N)
train <- df[1:split_idx, , drop = FALSE]
test  <- df[(split_idx + 1):N, , drop = FALSE]

train_X <- train %>% transmute(f1 = abs_ret_z, f2 = log_vol_z, f3 = park_21d_z, f4 = beta_60d_z)
test_X  <- test  %>% transmute(f1 = abs_ret_z, f2 = log_vol_z, f3 = park_21d_z, f4 = beta_60d_z)

#  Load saved fit/mapping OR fit on train 
if (USE_SAVED_ARTIFACTS && file.exists("hmm_fit.rds") && file.exists("state_map.csv")) {
  message("Loading saved HMM artifacts from Script 1…")
  fm <- readRDS("hmm_fit.rds")
  state_map_tbl <- readr::read_csv("state_map.csv", show_col_types = FALSE)
  state_map <- setNames(state_map_tbl$regime, as.character(state_map_tbl$state))
} else {
  message("Artifacts not found or disabled. Fitting on train (85%) and mapping states by train mean forward RV)…")
  set.seed(SEED)
  mod <- depmixS4::depmix(
    response = list(f1 ~ 1, f2 ~ 1, f3 ~ 1, f4 ~ 1),
    data     = as.data.frame(train_X),
    nstates  = NSTATES,
    family   = list(gaussian(), gaussian(), gaussian(), gaussian())
  )
  fm <- depmixS4::fit(mod, emcontrol = depmixS4::em.control(maxit = MAXIT), verbose = FALSE)
  post_train <- depmixS4::posterior(fm, type = "smoothing")
  W <- get_state_prob_matrix(post_train, NSTATES)
  rv_idx <- which(!is.na(train$rv21_fwd))
  state_means <- as.numeric(colSums(W[rv_idx,,drop=FALSE] * train$rv21_fwd[rv_idx]) / colSums(W[rv_idx,,drop=FALSE]))
  ord_states <- order(state_means)                                        # ascending mean RV
  state_map <- setNames(c("Low","Medium","High"), as.character(ord_states))
}

# Decode TEST with frozen parameters
post_test <- try_posterior(fm, as.data.frame(test_X))
prob_state_test <- get_state_prob_matrix(post_test, NSTATES)

# Map state probabilites -> regime probabilites using state_map order
ord <- as.integer(names(sort(setNames(1:3, names(state_map))[names(state_map)]))) # indices in Low/Med/High order
# Build permutation to L/M/H
perm <- match(c("Low","Medium","High"), state_map[as.character(1:NSTATES)])
prob_regime_test <- prob_state_test[, perm, drop = FALSE]
colnames(prob_regime_test) <- c("Low","Medium","High")

# Hard states & smoothing identical to Script 1
states_test <- try_viterbi(fm, as.data.frame(test_X))
states_test_smooth <- enforce_min_dwell(states_test, STICKY_DWELL)
regime_test <- factor(state_map[as.character(states_test_smooth)],
                      levels = c("Low","Medium","High"))

# Continuous risk score: E[regime index]
score <- prob_regime_test %*% matrix(c(1,2,3), ncol = 1)
score <- as.numeric(score)

#  Continuous metrics on last 15% 
ok <- which(!is.na(test$rv21_fwd))
y  <- test$rv21_fwd[ok]
s  <- score[ok]
pH <- prob_regime_test[ok, "High"]

rho <- suppressWarnings(cor(s, y, method = "spearman"))
tau <- suppressWarnings(cor(s, y, method = "kendall"))
r2  <- summary(lm(y ~ s))$r.squared

# Lift by score deciles
lift_tbl <- tibble(score = s, rv = y) %>%
  mutate(decile = ntile(score, 10)) %>%
  group_by(decile) %>%
  summarise(n = n(),
            mean_score = mean(score),
            mean_rv = mean(rv),
            median_rv = median(rv),
            .groups = "drop") %>%
  arrange(decile)

# Top-q capture (no binning of target)
q <- 0.10
thr <- quantile(s, probs = 1 - q, na.rm = TRUE)
is_top <- s >= thr
mean_top  <- mean(y[is_top])
mean_rest <- mean(y[!is_top])
uplift    <- mean_top - mean_rest
wt <- suppressWarnings(wilcox.test(y[is_top], y[!is_top], alternative = "greater"))

# Calibration by P(High) deciles
calib_tbl <- tibble(pHigh = as.numeric(pH), rv = y) %>%
  mutate(bin = ntile(pHigh, 10)) %>%
  group_by(bin) %>%
  summarise(n = n(), mean_pHigh = mean(pHigh), mean_rv = mean(rv),
            median_rv = median(rv), .groups = "drop") %>%
  arrange(bin)

q_hi <- 0.20
hi_cut <- quantile(train$rv21_fwd, 1 - q_hi, na.rm = TRUE)   # train-only threshold
y_true_bin <- as.integer(test$rv21_fwd[ok] >= hi_cut)        # 1 = High, 0 = Rest
p_high <- as.numeric(prob_regime_test[ok, "High"])


cat("==== OOS metrics on last 15% ====\n")
cat(sprintf("Spearman rho(score, future RV21): %.4f\n", rho))
cat(sprintf("ROC-AUC (High vs Rest, top %.0f%% by train RV): %.4f\n", 100*q_hi, as.numeric(auc_val)))

# Save
df$set <- c(rep("train", nrow(train)), rep("test", nrow(test)))
readr::write_csv(df, OUT_TIMESERIES_CSV)

oos_out <- tibble(
  date   = test$date,
  regime = regime_test,
  state  = as.integer(states_test_smooth),
  score  = score
) %>%
  bind_cols(as_tibble(prob_regime_test)) %>%
  rename(prob_Low = Low, prob_Medium = Medium, prob_High = High)

readr::write_csv(oos_out, OUT_PRED_CSV)
readr::write_csv(lift_tbl, OUT_LIFT_CSV)
readr::write_csv(calib_tbl, OUT_CALIB_CSV)

metrics_txt <- sprintf(paste(
  "Continuous OOS metrics on last %.0f%% (no tertiles)\n",
  "Spearman rho(score, RV21_fwd) = %.4f\n",
  "Kendall tau(score, RV21_fwd)  = %.4f\n",
  "Linear R^2 (RV21_fwd ~ score) = %.4f\n",
  "Top %.0f%% score capture: mean_top=%.6f | mean_rest=%.6f | uplift=%.6f | Wilcoxon p=%.3g\n",
  sep = ""),
  100*OOS_FRAC, rho, tau, r2, 100*q, mean_top, mean_rest, uplift, wt$p.value
)
writeLines(metrics_txt, OUT_METRICS_TXT)

message(sprintf("Saved: %s, %s, %s, %s, %s",
                OUT_TIMESERIES_CSV, OUT_PRED_CSV, OUT_LIFT_CSV, OUT_CALIB_CSV, OUT_METRICS_TXT))

#Spit Out Spearman Figure ISOLATED
# Spearman rho between score and forward realized vol
spearman_rho <- suppressWarnings(cor(s, y, method = "spearman", use = "complete.obs"))
cat(sprintf("Spearman rho(score, future RV21) = %.4f\n", spearman_rho))

} # End of if(packages_loaded) block

# Export posterior probabilities P(z|t) for plotting - MINIMAL ADDITION
# Try to use real probabilities if available, otherwise use synthetic
if(exists("prob_regime_test") && exists("test")) {
  rgb_plot_data <- data.frame(
    date = test$date,
    price = test$adjusted,
    prob_R = prob_regime_test[, "High"],     # High vol = Red
    prob_G = prob_regime_test[, "Low"],      # Low vol = Green
    prob_B = prob_regime_test[, "Medium"]    # Medium vol = Blue
  )
  write.csv(rgb_plot_data, "hmm_rgb_data.csv", row.names = FALSE)
  message("Exported P(z|t) to hmm_rgb_data.csv for RGB plotting")
} else {
  # SYNTHETIC FALLBACK - Create demo data if packages failed
  cat("Creating synthetic RGB data for demo...\n")
  n_days <- 200
  dates <- seq(as.Date("2022-01-01"), length.out = n_days, by = "day")
  prices <- 100 * cumprod(1 + rnorm(n_days, 0, 0.02))

  # Generate synthetic probabilities with regime persistence
  prob_R <- prob_G <- prob_B <- numeric(n_days)
  current_state <- sample(1:3, 1)
  for(i in 1:n_days) {
    if(runif(1) < 0.05) current_state <- sample(1:3, 1)  # 5% chance of switch
    probs <- c(0.1, 0.1, 0.1)
    probs[current_state] <- 0.8
    prob_R[i] <- probs[1]; prob_G[i] <- probs[2]; prob_B[i] <- probs[3]
  }

  write.csv(data.frame(date=dates, price=prices, prob_R=prob_R, prob_G=prob_G, prob_B=prob_B),
            "hmm_rgb_data.csv", row.names=FALSE)
  message("Created synthetic hmm_rgb_data.csv for RGB plotting demo")
}

