############################################################
# TABLE 6-STYLE PANDEMIC BOND PRICING
# SIRD replacement of the Risks Table 6 model
#
# New trigger:
#   D_total(t) >= 2500
#   maI(t)     >= 250
#
# where
#   C_total(t)            = I(t) + R(t) + D(t)
#   daily_new_infections  = max(0, C_total(t) - C_total(t-1))
#   maI(t)                = 7-day moving average of daily new infections
#
# This script keeps the pricing shell from Risks:
#   - 36 predicted coupons from Table 3
#   - face value = 225M
#   - required yield = 8.6734%
#   - P(H) = 0.1393
#   - M = 5000 simulation paths
#   - volatility sweep 0.04, 0.08, ..., 0.40
#
# but replaces the epidemic model + trigger by SIRD.
############################################################

#############################
# 0) GLOBAL SETTINGS
#############################

options(stringsAsFactors = FALSE)

#############################
# 1) BOND INPUTS FROM RISKS
#############################

BOND <- list(
  face_value = 225e6,
  y_required = 0.086734,   # 8.6734%
  P_H        = 0.1393,     # probability of pandemic
  M_paths    = 5000
)

# Coupon payment dates from Table 3
coupon_dates <- as.Date(c(
  "2017-08-15","2017-09-15","2017-10-16","2017-11-15","2017-12-15","2018-01-15",
  "2018-02-15","2018-03-15","2018-04-16","2018-05-15","2018-06-15","2018-07-16",
  "2018-08-15","2018-09-17","2018-10-15","2018-11-15","2018-12-17","2019-01-15",
  "2019-02-15","2019-03-15","2019-04-15","2019-05-15","2019-06-17","2019-07-15",
  "2019-08-15","2019-09-16","2019-10-15","2019-11-15","2019-12-16","2020-01-15",
  "2020-02-17","2020-03-16","2020-04-15","2020-05-15","2020-06-15","2020-07-15"
))

issue_date    <- as.Date("2017-07-07")
maturity_date <- as.Date("2020-07-15")

# Days from issue date to coupon dates
coupon_days <- as.integer(coupon_dates - issue_date) + 1L

# Predicted coupon amounts from Table 3, in USD millions
coupon_amounts_million <- c(
  1.941576, 1.552743, 1.558080, 1.508914, 1.506405, 1.550377,
  1.546438, 1.403378, 1.611381, 1.467906, 1.578030, 1.588832,
  1.547175, 1.708552, 1.455233, 1.616731, 1.675545, 1.524186,
  1.633633, 1.478876, 1.640417, 1.590558, 1.752495, 1.489606,
  1.652194, 1.710977, 1.555140, 1.667228, 1.672601, 1.623373,
  1.789533, 1.520603, 1.631066, 1.632746, 1.688736, 1.635570
)

coupon_amounts <- coupon_amounts_million * 1e6

#############################
# 2) SIRD MODEL CONFIG
#############################

CFG_BASE <- list(
  # horizon / discretization
  T_end = max(coupon_days),   # 1104 days
  dt    = 1,                  # daily step, Table 6-style horizon
  M     = BOND$M_paths,       # 5000 paths
  
  # population
  N_pop = 8.30e9,
  
  # initial counts
  I0_count = 100,
  R0_count = 0,
  D0_count = 0,
  
  # SIRD-OU parameters
  beta  = 0.22,
  beta0 = 0.15,
  theta = 0.13,
  xi    = 0.12,   # this will be varied in the Table-6 style volatility sweep
  alpha = 0.10,
  gamma = 0.09,
  
  # trigger thresholds
  theta_D_total = 2500,
  theta_I_ma7   = 250,
  
  # reproducibility
  seed = 123
)

#############################
# 3) HELPERS
#############################

beta_bar <- function(t, beta, beta0, theta) {
  beta + (beta0 - beta) * exp(-theta * t)
}

sigma_t <- function(t, theta, xi) {
  (xi / sqrt(2 * theta)) * sqrt(1 - exp(-2 * theta * t))
}

roll_mean_7_matrix <- function(X) {
  # X: n x M matrix
  n <- nrow(X)
  M <- ncol(X)
  out <- matrix(NA_real_, nrow = n, ncol = M)
  if (n < 7) return(out)
  
  for (i in 7:n) {
    out[i, ] <- colMeans(X[(i - 6):i, , drop = FALSE])
  }
  out
}

discount_cashflows <- function(amounts, days, y) {
  sum(amounts / (1 + y)^(days / 360))
}

bond_value_full <- function(bond, coupon_days, coupon_amounts) {
  discount_cashflows(coupon_amounts, coupon_days, bond$y_required) +
    bond$face_value / (1 + bond$y_required)^(max(coupon_days) / 360)
}

bond_value_triggered <- function(tau_day, bond, coupon_days, coupon_amounts) {
  # Coupons strictly before trigger day are received.
  keep <- which(coupon_days < tau_day)
  
  if (length(keep) == 0) {
    return(0)
  }
  
  discount_cashflows(coupon_amounts[keep], coupon_days[keep], bond$y_required)
}

#############################
# 4) SIMULATE DAILY SIRD PATHS
#
# This is a daily-step version intended for the Table 6
# experiment with 5000 paths and 1104 days.
#############################

simulate_sird_paths_daily <- function(cfg) {
  
  if (!is.null(cfg$seed)) set.seed(cfg$seed)
  
  T_end <- cfg$T_end
  M     <- cfg$M
  N_pop <- cfg$N_pop
  
  beta  <- cfg$beta
  beta0 <- cfg$beta0
  theta <- cfg$theta
  xi    <- cfg$xi
  alpha <- cfg$alpha
  gamma <- cfg$gamma
  
  # Initial fractions
  I0 <- cfg$I0_count / N_pop
  R0 <- cfg$R0_count / N_pop
  D0 <- cfg$D0_count / N_pop
  S0 <- 1 - I0 - R0 - D0
  
  if (S0 <= 0) stop("Initial S0 must be positive.")
  
  days <- 0:T_end
  n    <- length(days)
  
  S <- matrix(NA_real_, nrow = n, ncol = M)
  I <- matrix(NA_real_, nrow = n, ncol = M)
  R <- matrix(NA_real_, nrow = n, ncol = M)
  D <- matrix(NA_real_, nrow = n, ncol = M)
  
  S[1, ] <- S0
  I[1, ] <- I0
  R[1, ] <- R0
  D[1, ] <- D0
  
  for (k in 1:T_end) {
    
    tk  <- days[k]
    btk <- beta_bar(tk, beta, beta0, theta)
    sig <- sigma_t(tk, theta, xi)
    
    Zk  <- rnorm(M)
    dBk <- Zk   # because dt = 1, sqrt(dt) = 1
    
    S_prev <- S[k, ]
    I_prev <- I[k, ]
    R_prev <- R[k, ]
    D_prev <- D[k, ]
    
    # stochastic log-Euler for I and S over one day
    eta_I <- sig * S_prev
    drift_logI <- btk * S_prev - (alpha + gamma) - 0.5 * eta_I^2
    I_new <- I_prev * exp(drift_logI + eta_I * dBk)
    
    eta_S <- -sig * I_prev
    drift_logS <- -btk * I_prev - 0.5 * eta_S^2
    S_new <- S_prev * exp(drift_logS + eta_S * dBk)
    
    # Euler updates for R and D
    R_new <- R_prev + alpha * I_prev
    D_new <- D_prev + gamma * I_prev
    
    # Non-negativity cleanup
    S_new <- pmax(S_new, 0)
    I_new <- pmax(I_new, 0)
    R_new <- pmax(R_new, 0)
    D_new <- pmax(D_new, 0)
    
    # Optional mass normalization if total creeps above 1
    total_new <- S_new + I_new + R_new + D_new
    bad <- total_new > 1
    if (any(bad)) {
      S_new[bad] <- S_new[bad] / total_new[bad]
      I_new[bad] <- I_new[bad] / total_new[bad]
      R_new[bad] <- R_new[bad] / total_new[bad]
      D_new[bad] <- D_new[bad] / total_new[bad]
    }
    
    S[k + 1, ] <- S_new
    I[k + 1, ] <- I_new
    R[k + 1, ] <- R_new
    D[k + 1, ] <- D_new
  }
  
  list(
    days = days,
    S = S,
    I = I,
    R = R,
    D = D
  )
}

#############################
# 5) BUILD SIRD TRIGGER SERIES
#############################

build_trigger_series_sird <- function(sim, cfg) {
  
  N_pop <- cfg$N_pop
  
  I_count <- sim$I * N_pop
  R_count <- sim$R * N_pop
  D_total <- sim$D * N_pop
  
  # cumulative ever-infected
  C_total <- I_count + R_count + D_total
  
  n <- nrow(C_total)
  M <- ncol(C_total)
  
  daily_new_inf <- matrix(NA_real_, nrow = n, ncol = M)
  daily_new_inf[1, ] <- NA_real_
  
  if (n >= 2) {
    daily_new_inf[2:n, ] <- pmax(0, C_total[2:n, , drop = FALSE] - C_total[1:(n - 1), , drop = FALSE])
  }
  
  maI <- roll_mean_7_matrix(daily_new_inf)
  
  list(
    I_count       = I_count,
    R_count       = R_count,
    D_total       = D_total,
    C_total       = C_total,
    daily_new_inf = daily_new_inf,
    maI           = maI
  )
}

#############################
# 6) CALCULATE TRIGGER TIME
#############################

calc_trigger_tau_sird <- function(days, D_total, maI, cfg) {
  
  M <- ncol(D_total)
  tau <- rep(NA_integer_, M)
  
  for (j in 1:M) {
    cond <- (!is.na(D_total[, j])) &
      (!is.na(maI[, j])) &
      (D_total[, j] >= cfg$theta_D_total) &
      (maI[, j]     >= cfg$theta_I_ma7)
    
    hit <- which(cond)
    if (length(hit) > 0) {
      tau[j] <- days[hit[1]]
    }
  }
  
  tau
}

#############################
# 7) RUN ONE VOLATILITY LEVEL
#############################

run_one_xi <- function(xi_value,
                       cfg_base,
                       bond,
                       coupon_days,
                       coupon_amounts) {
  
  cfg <- cfg_base
  cfg$xi <- xi_value
  
  sim <- simulate_sird_paths_daily(cfg)
  series <- build_trigger_series_sird(sim, cfg)
  
  tau <- calc_trigger_tau_sird(
    days    = sim$days,
    D_total = series$D_total,
    maI     = series$maI,
    cfg     = cfg
  )
  
  triggered <- is.finite(tau)
  
  # Conditional trigger probability from simulation
  p_cond <- mean(triggered)
  
  # Unconditional trigger probability, same shell as Risks
  P_C <- p_cond * bond$P_H
  
  BV_full <- bond_value_full(bond, coupon_days, coupon_amounts)
  
  if (any(triggered)) {
    BV_partial <- vapply(
      tau[triggered],
      FUN = bond_value_triggered,
      FUN.VALUE = numeric(1),
      bond = bond,
      coupon_days = coupon_days,
      coupon_amounts = coupon_amounts
    )
    E_trigger <- mean(BV_partial)
  } else {
    E_trigger <- BV_full
  }
  
  BondPrice <- E_trigger * P_C + BV_full * (1 - P_C)
  
  data.frame(
    xi           = xi_value,
    Bond.Price   = BondPrice / 1e6,  # report in USD millions
    Prob.Trigger = P_C
  )
}

#############################
# 8) BUILD TABLE 6-STYLE OUTPUT
#############################

make_table6_sird <- function(cfg_base,
                             bond,
                             coupon_days,
                             coupon_amounts,
                             xi_grid = seq(0.04, 0.40, by = 0.04)) {
  
  out <- do.call(
    rbind,
    lapply(
      xi_grid,
      run_one_xi,
      cfg_base = cfg_base,
      bond = bond,
      coupon_days = coupon_days,
      coupon_amounts = coupon_amounts
    )
  )
  
  rownames(out) <- NULL
  out
}

#############################
# 9) OPTIONAL SINGLE-PATH PLOT
#############################

plot_trigger_all_in_one_sird <- function(days, D_total, maI, tau, cfg) {
  
  M <- ncol(D_total)
  col_D   <- "orange"
  col_maI <- "blue"
  
  yD_rng   <- range(D_total, na.rm = TRUE)
  ymaI_rng <- range(maI, na.rm = TRUE)
  
  if (!all(is.finite(yD_rng)))   yD_rng   <- c(0, 1)
  if (!all(is.finite(ymaI_rng))) ymaI_rng <- c(0, 1)
  
  if (diff(yD_rng) == 0)   yD_rng   <- yD_rng + c(-0.5, 0.5)
  if (diff(ymaI_rng) == 0) ymaI_rng <- ymaI_rng + c(-0.5, 0.5)
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mar = c(5, 5, 4, 6))
  
  # left axis
  matplot(days, D_total, type = "l", lty = 1, lwd = 2,
          col = col_D,
          xlab = "Days",
          ylab = "Total deaths D_total(t) [left axis]",
          main = paste0("SIRD Trigger Graph: D_total (orange), maI (blue), M = ", M),
          ylim = yD_rng)
  
  abline(h = cfg$theta_D_total, col = col_D, lty = 3, lwd = 2)
  
  tau_ok <- tau[is.finite(tau)]
  if (length(tau_ok) > 0) {
    abline(v = tau_ok, col = "red", lty = 2, lwd = 1)
  }
  
  for (j in 1:M) {
    if (is.finite(tau[j])) {
      idx <- which(days == tau[j])
      if (length(idx) == 1) {
        points(days[idx], D_total[idx, j], pch = 17, cex = 1.2, col = "darkred")
      }
    }
  }
  
  # right axis
  par(new = TRUE)
  matplot(days, maI, type = "l", lty = 1, lwd = 2,
          col = col_maI,
          axes = FALSE, xlab = "", ylab = "",
          ylim = ymaI_rng)
  
  axis(side = 4, col = col_maI, col.axis = col_maI)
  mtext("maI(t): 7-day moving average of daily new infections [right axis]",
        side = 4, line = 2, col = col_maI)
  
  abline(h = cfg$theta_I_ma7, col = col_maI, lty = 3, lwd = 2)
  
  for (j in 1:M) {
    if (is.finite(tau[j])) {
      idx <- which(days == tau[j])
      if (length(idx) == 1) {
        points(days[idx], maI[idx, j], pch = 19, cex = 1.0, col = "red")
      }
    }
  }
  
  legend("topleft",
         legend = c("D_total paths", "maI paths", "death threshold = 2500",
                    "maI threshold = 250", "trigger tau",
                    "trigger point on maI", "trigger point on D_total"),
         col    = c(col_D, col_maI, col_D, col_maI, "red", "red", "darkred"),
         lty    = c(1, 1, 3, 3, 2, NA, NA),
         lwd    = c(2, 2, 2, 2, 1, NA, NA),
         pch    = c(NA, NA, NA, NA, NA, 19, 17),
         bty    = "n")
}

#############################
# 10) RUN THE TABLE 6-STYLE SIRD STUDY
#############################

xi_grid <- seq(0.04, 0.40, by = 0.04)

table6_sird <- make_table6_sird(
  cfg_base       = CFG_BASE,
  bond           = BOND,
  coupon_days    = coupon_days,
  coupon_amounts = coupon_amounts,
  xi_grid        = xi_grid
)

# Pretty print
table6_sird_print <- table6_sird
table6_sird_print$Bond.Price   <- round(table6_sird_print$Bond.Price, 2)
table6_sird_print$Prob.Trigger <- round(table6_sird_print$Prob.Trigger, 4)

print(table6_sird_print)

#############################
# 11) OPTIONAL: ONE SAMPLE PATH CHECK
#############################

# Smaller run for plotting / debugging
CFG_CHECK <- CFG_BASE
CFG_CHECK$M    <- 1
CFG_CHECK$xi   <- 0.12
CFG_CHECK$seed <- 123

sim_check    <- simulate_sird_paths_daily(CFG_CHECK)
series_check <- build_trigger_series_sird(sim_check, CFG_CHECK)
tau_check    <- calc_trigger_tau_sird(sim_check$days,
                                      series_check$D_total,
                                      series_check$maI,
                                      CFG_CHECK)

cat("\n========================================\n")
cat("SIRD CHECK RUN\n")
cat("========================================\n")
cat("tau =", tau_check, "\n")
if (is.finite(tau_check[1])) {
  idx <- which(sim_check$days == tau_check[1])
  cat("D_total(tau) =", series_check$D_total[idx, 1], "\n")
  cat("maI(tau)     =", series_check$maI[idx, 1], "\n")
} else {
  cat("No trigger before maturity.\n")
}

plot_trigger_all_in_one_sird(
  days    = sim_check$days,
  D_total = series_check$D_total,
  maI     = series_check$maI,
  tau     = tau_check,
  cfg     = CFG_CHECK
)



# ============================================================
# FORMAT THE SIRD RESULTS LIKE TABLE 6 IN THE PAPER
# ============================================================


# table6_sird is assumed to come from:
# table6_sird <- make_table6_sird(...)

# extract vectors in the same style as your old code
xi_target  <- table6_sird$xi
Bond_Price <- table6_sird$Bond.Price
P_C        <- table6_sird$Prob.Trigger

# optional ordering, in case you want to be explicit
idx <- order(xi_target)

tab6 <- rbind(
  Bond.Price   = round(Bond_Price[idx], 2),
  Prob.Trigger = round(P_C[idx], 2)
)
# 2-row matrix: bond prices and trigger probabilities

colnames(tab6) <- sprintf("%.2f", xi_target[idx])
# column names: "0.04", ..., "0.40"

knitr::kable(
  tab6,
  rownames = TRUE,
  col.names = c("xi", sprintf("%.2f", xi_target[idx])),
  align = c("l", rep("c", ncol(tab6))),
  caption = "Table 6. The bond price movement against volatility terms."
)
# Prints a table with the top-left header cell as "xi"

mean(Bond_Price)
# average of all Bond_Price values across xi grid (not in the paper table)


# ---------------------------
# 7) Save outputs
# ---------------------------

saveRDS(tab6, file = "table6_sird_like_paper.rds")
write.csv(tab6, file = "table6_sird_like_paper.csv")