#Title: r2k_hmm_project
suppressPackageStartupMessages({
  pkgs <- c("tidyverse","lubridate","tidyquant","zoo","slider",
            "PerformanceAnalytics","depmixS4")
  to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
  if(length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
  lapply(pkgs, library, character.only = TRUE)
})

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
  # returns a vector same length as x/y; NA for first n-1 entries(?)
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
  # Try depmixS4::viterbi if available; otherwise return MAP from posterior()
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

# Optional: light hysteresis to cut thrash (min dwell time)
apply_hysteresis <- function(states, min_dwell = 3) {
  if (all(is.na(states))) return(states)
  s <- states
  r <- rle(s)
  r$lengths[r$lengths < min_dwell] <- min_dwell
  inverse.rle(r)[seq_along(s)]
}

# Data download, scraping, and bucketing done below: (Ask AJay for Help)
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

# Train / Validate:
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

# Fit HMM 
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
# Use fitted params on the full data to get a single, consistent state mapping
states_raw <- try_viterbi(fm, as.data.frame(all_X))

# Optional: apply gentle hysteresis to reduce day-to-day churn
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


