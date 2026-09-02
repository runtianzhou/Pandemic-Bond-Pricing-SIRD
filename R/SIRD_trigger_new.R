############################################################
# SIRD-OU simulation + custom trigger
#
# Trigger:
#   maI(t)     >= 250
#   D_total(t) >= 2500
#
# IMPORTANT:
# Here maI(t) means:
#   7-day moving average of DAILY NEW INFECTIONS
# not the moving average of active infected I(t).
#
# Definitions used in this code:
#   cumulative_cases(t)     = I(t) + R(t) + D(t)
#   daily_new_infections(t) = max(0, cumulative_cases(t) - cumulative_cases(t-1))
#   maI(t)                  = 7-day moving average of daily_new_infections
#   D_total(t)              = cumulative deaths
############################################################

#############################
# 0) CONFIG
#############################
CFG <- list(
  # Total simulation horizon in days
  T_end = 200,
  
  # Output step in days
  # Here dt = 1 means we keep daily observations
  dt    = 1,
  
  # Fine simulation step inside each day
  # Since dt_n = 0.01, each day is broken into 100 fine steps
  dt_n  = 1,
  
  # Number of simulated paths
  M     = 1,
  
  # Total population size used to convert fractions to counts
  N_pop = 8.30e9,
  
  # Initial counts at time 0
  I0_count = 100,
  R0_count = 0,
  D0_count = 0,
  
  # SIRD-OU model parameters
  beta  = 0.22,
  beta0 = 0.15,
  theta = 0.13,
  xi    = 0.12,
  alpha = 0.10,
  gamma = 0.09,
  
  # Trigger thresholds:
  # cumulative deaths must be at least 2500
  theta_D_total = 2500,
  
  # 7-day moving average of DAILY NEW infections must be at least 250
  theta_I_ma7   = 250,
  
  # Plot window
  x_min = 0,
  x_max = 200,
  
  # Random seed; NULL means no fixed seed
  seed = NULL
)

#############################
# 1) HELPERS
#############################

# Time-varying mean-reverting beta-bar function
# This matches the deterministic part of the time-varying transmission structure
beta_bar <- function(t, beta, beta0, theta) {
  beta + (beta0 - beta) * exp(-theta * t)
}

# Time-varying volatility scale used in the stochastic part
sigma_t <- function(t, theta, xi) {
  (xi / sqrt(2 * theta)) * sqrt(1 - exp(-2 * theta * t))
}

# Checks that dt_n exactly partitions dt
# Example: dt = 1 and dt_n = 0.01 gives m = 100 fine steps per day
check_partition <- function(dt, dt_n) {
  m_real <- dt / dt_n
  m <- as.integer(round(m_real))
  if (abs(m_real - m) > 1e-12) {
    stop("dt_n must partition dt exactly.")
  }
  m
}

# Computes a trailing 7-day moving average
# For day i, it averages x[i-6], ..., x[i]
# First 6 positions are NA because there are not yet 7 days available
roll_mean_7 <- function(x) {
  n <- length(x)
  out <- rep(NA_real_, n)
  if (n < 7) return(out)
  for (i in 7:n) {
    out[i] <- mean(x[(i - 6):i], na.rm = FALSE)
  }
  out
}

#############################
# 2) SIMULATE SIRD PATHS
#############################
simulate_sird_paths <- function(cfg) {
  # Read configuration values into local variables
  T_end <- cfg$T_end
  dt    <- cfg$dt
  dt_n  <- cfg$dt_n
  M     <- cfg$M
  N_pop <- cfg$N_pop
  
  beta  <- cfg$beta
  beta0 <- cfg$beta0
  theta <- cfg$theta
  xi    <- cfg$xi
  alpha <- cfg$alpha
  gamma <- cfg$gamma
  
  # Set random seed if provided
  if (!is.null(cfg$seed)) set.seed(cfg$seed)
  
  # Convert initial counts to fractions of total population
  I0 <- cfg$I0_count / N_pop
  R0 <- cfg$R0_count / N_pop
  D0 <- cfg$D0_count / N_pop
  
  # Susceptible fraction is whatever is left
  S0 <- 1 - I0 - R0 - D0
  
  # Sanity check: susceptible fraction must be positive
  if (S0 <= 0) stop("Initial S0 must be positive.")
  
  # Number of fine steps per day
  m <- check_partition(dt, dt_n)
  
  # Number of whole days in the horizon
  n_days_real <- T_end / dt
  n_days <- as.integer(round(n_days_real))
  
  # Sanity check: horizon must be an integer number of days
  if (abs(n_days_real - n_days) > 1e-12) {
    stop("T_end must be an integer number of days when dt = 1.")
  }
  
  # Total number of fine simulation steps
  n_steps <- n_days * m
  
  # Fine time grid: 0, dt_n, 2*dt_n, ...
  t_fine  <- seq(from = 0, by = dt_n, length.out = n_steps + 1)
  
  # Allocate matrices to store fine-grid paths
  S <- matrix(NA_real_, nrow = n_steps + 1, ncol = M)
  I <- matrix(NA_real_, nrow = n_steps + 1, ncol = M)
  R <- matrix(NA_real_, nrow = n_steps + 1, ncol = M)
  D <- matrix(NA_real_, nrow = n_steps + 1, ncol = M)
  
  # Set initial values at time 0
  S[1, ] <- S0
  I[1, ] <- I0
  R[1, ] <- R0
  D[1, ] <- D0
  
  # Fine-step simulation loop
  for (k in 1:n_steps) {
    # Current fine-grid time
    tk  <- t_fine[k]
    
    # Time-varying transmission drift component
    btk <- beta_bar(tk, beta, beta0, theta)
    
    # Time-varying volatility
    sig <- sigma_t(tk, theta, xi)
    
    # One standard normal shock per path for this fine step
    Zk  <- rnorm(M)
    
    # Brownian increment over length dt_n
    dBk <- sqrt(dt_n) * Zk
    
    # Previous state values
    S_prev <- S[k, ]
    I_prev <- I[k, ]
    R_prev <- R[k, ]
    D_prev <- D[k, ]
    
    # ---- Infectious process I ----
    eta_I <- sig * S_prev
    drift_logI <- btk * S_prev - (alpha + gamma) - 0.5 * (eta_I^2)
    I_new <- I_prev * exp(drift_logI * dt_n + eta_I * dBk)
    
    # ---- Susceptible process S ----
    eta_S <- -sig * I_prev
    drift_logS <- -btk * I_prev - 0.5 * (eta_S^2)
    S_new <- S_prev * exp(drift_logS * dt_n + eta_S * dBk)
    
    # Numerical cleanup for S and I first
    S_new <- pmax(S_new, 0)
    I_new <- pmax(I_new, 0)
    
    # If needed, cap S + I at 1 to preserve fraction interpretation
    total_SI <- S_new + I_new
    over <- total_SI > 1
    if (any(over)) {
      S_new[over] <- S_new[over] / total_SI[over]
      I_new[over] <- I_new[over] / total_SI[over]
    }
    
    # ---- Recovered and dead processes R and D ----
    remaining_RD <- 1 - S_new - I_new
    R_new <- (alpha / (alpha + gamma)) * remaining_RD
    D_new <- (gamma / (alpha + gamma)) * remaining_RD
    
    # Numerical cleanup: prevent negative values
    S_new <- pmax(S_new, 0)
    I_new <- pmax(I_new, 0)
    R_new <- pmax(R_new, 0)
    D_new <- pmax(D_new, 0)
    
    # Store next-step values
    S[k + 1, ] <- S_new
    I[k + 1, ] <- I_new
    R[k + 1, ] <- R_new
    D[k + 1, ] <- D_new
  }
  
  # Keep only daily endpoints from the fine grid
  idx_day_end <- 1 + (0:n_days) * m
  days <- 0:n_days
  
  # Daily sampled paths
  S_day <- S[idx_day_end, , drop = FALSE]
  I_day <- I[idx_day_end, , drop = FALSE]
  R_day <- R[idx_day_end, , drop = FALSE]
  D_day <- D[idx_day_end, , drop = FALSE]
  
  # Return daily paths and initial fractions
  list(
    S = S_day,
    I = I_day,
    R = R_day,
    D = D_day,
    days = days,
    init_fraction = c(S0 = S0, I0 = I0, R0 = R0, D0 = D0)
  )
}

#############################
# 3) BUILD TRIGGER SERIES
#############################
build_trigger_series_sird <- function(sim, cfg) {
  # Daily time vector
  days  <- sim$days
  
  # Number of paths
  M     <- ncol(sim$I)
  
  # Population size for converting fractions to counts
  N_pop <- cfg$N_pop
  
  # Convert daily fractions to daily counts
  I_count <- sim$I * N_pop
  R_count <- sim$R * N_pop
  D_total <- sim$D * N_pop
  
  # IMPORTANT:
  # This defines cumulative ever-infected as I + R + D
  # That is correct in a SIRD model because everyone who has ever been infected
  # is either currently infected, recovered, or dead.
  C_total <- I_count + R_count + D_total
  
  # Number of daily observations
  n <- length(days)
  
  # Allocate containers
  daily_new_inf <- matrix(NA_real_, nrow = n, ncol = M)
  maI           <- matrix(NA_real_, nrow = n, ncol = M)
  
  # Build daily new infections and 7-day moving average path by path
  for (j in 1:M) {
    # Day 0 has no previous day, so daily increment is undefined
    daily_new_inf[1, j] <- NA_real_
    
    for (i in 2:n) {
      # Daily new infections = increase in cumulative ever-infected
      # max(0, ...) is used to prevent tiny negative numerical artifacts
      daily_new_inf[i, j] <- max(0, C_total[i, j] - C_total[i - 1, j])
    }
    
    # 7-day moving average of DAILY NEW INFECTIONS
    # This is exactly your intended maI trigger series
    maI[, j] <- roll_mean_7(daily_new_inf[, j])
  }
  
  # Return all trigger-related series
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
# 4) TRIGGER TIME PER PATH
#############################
calc_trigger_tau_sird <- function(days, D_total, maI, cfg) {
  # Number of paths
  M <- ncol(D_total)
  
  # tau[j] = first trigger day for path j
  tau <- rep(NA_integer_, M)
  
  for (j in 1:M) {
    # Trigger condition:
    # 1) deaths must be available and maI must be available
    # 2) cumulative deaths >= 2500
    # 3) 7-day MA of daily new infections >= 250
    cond <- (!is.na(D_total[, j])) & (!is.na(maI[, j])) &
      (D_total[, j] >= cfg$theta_D_total) &
      (maI[, j]     >= cfg$theta_I_ma7)
    
    # Find all days where both trigger conditions hold
    hit <- which(cond)
    
    # Store the FIRST such day
    if (length(hit) > 0) tau[j] <- days[hit[1]]
  }
  
  tau
}

#############################
# 5) PLOT ALL IN ONE
#    left axis  = D_total
#    right axis = maI
#############################
plot_trigger_all_in_one_sird <- function(days, D_total, maI, tau, cfg) {
  x_min <- cfg$x_min
  x_max <- cfg$x_max
  M     <- ncol(D_total)
  
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
  
  # ----------------------------
  # LEFT AXIS: D_total
  # ----------------------------
  matplot(days, D_total, type = "l", lty = 1, lwd = 2,
          col = col_D,
          xlab = "Days",
          ylab = "Total confirmed deaths D_total(t) [left axis]",
          main = paste0("SIRD Trigger Graph: D_total (orange), maI (blue), M = ", M),
          xlim = c(x_min, x_max),
          ylim = yD_rng)
  
  abline(h = cfg$theta_D_total, col = col_D, lty = 4, lwd = 2)
  
  tau_ok <- tau[is.finite(tau)]
  if (length(tau_ok) > 0) {
    abline(v = tau_ok, col = "red", lty = 4, lwd = 2)
  }
  
  # ORANGE trigger point must be added here, on left-axis layer
  for (j in 1:M) {
    if (is.finite(tau[j])) {
      idx <- which(days == tau[j])
      if (length(idx) == 1) {
        points(days[idx], D_total[idx, j],
               pch = 17, cex = 1.3, col = "darkred")
      }
    }
  }
  
  # ----------------------------
  # RIGHT AXIS: maI
  # ----------------------------
  par(new = TRUE)
  matplot(days, maI, type = "l", lty = 1, lwd = 2,
          col = col_maI,
          axes = FALSE, xlab = "", ylab = "",
          xlim = c(x_min, x_max),
          ylim = ymaI_rng)
  
  axis(side = 4, col = col_maI, col.axis = col_maI)
  mtext("maI(t): 7-day moving average of daily new infections [right axis]",
        side = 4, line = 2, col = col_maI)
  
  abline(h = cfg$theta_I_ma7, col = col_maI, lty =4, lwd = 2)
  
  # RED trigger point on blue maI curve
  for (j in 1:M) {
    if (is.finite(tau[j])) {
      idx <- which(days == tau[j])
      if (length(idx) == 1) {
        points(days[idx], maI[idx, j],
               pch = 19, cex = 1.0, col = "red")
      }
    }
  }
  
  legend("topleft",
         legend = c("D_total paths", "maI paths", "death threshold = 2500",
                    "maI threshold = 250", "trigger tau",
                    "trigger point on maI", "trigger point on D_total"),
         col    = c(col_D, col_maI, col_D, col_maI, "red", "red", "darkred"),
         lty    = c(1, 1, 4, 4, 4, NA, NA),
         lwd    = c(2, 2, 2, 2, 2, NA, NA),
         pch    = c(NA, NA, NA, NA, NA, 19, 17),
         bty    = "n")
}

#############################
# 6) PRINT SUMMARY
#############################
print_trigger_summary_sird <- function(days, D_total, maI, tau, cfg) {
  M <- ncol(D_total)
  
  cat("========================================\n")
  cat("SIRD TRIGGER SUMMARY\n")
  cat("========================================\n")
  cat("Population N_pop             =", format(cfg$N_pop, scientific = TRUE), "\n")
  cat("Initial infected count       =", cfg$I0_count, "\n")
  cat("Death threshold              =", cfg$theta_D_total, "\n")
  cat("maI threshold                =", cfg$theta_I_ma7, "\n")
  cat("maI = 7-day avg daily new infections\n\n")
  
  cat("tau per path:\n")
  print(tau)
  cat("\n")
  
  for (j in 1:M) {
    cat("----------------------------------------\n")
    cat("Path", j, "\n")
    cat("----------------------------------------\n")
    
    if (is.finite(tau[j])) {
      idx <- which(days == tau[j])
      cat("Triggered at day tau =", tau[j], "\n")
      cat("D_total(tau)         =", D_total[idx, j], "\n")
      cat("maI(tau)             =", maI[idx, j], "\n\n")
    } else {
      cat("No trigger before horizon.\n\n")
    }
  }
}

#############################
# 7) RUN
#############################

# Simulate daily SIRD paths
sim <- simulate_sird_paths(CFG)

# Build the trigger-related series from the simulated paths
series  <- build_trigger_series_sird(sim, CFG)

# Extract cumulative deaths
D_total <- series$D_total

# Extract 7-day moving average of DAILY NEW INFECTIONS
maI     <- series$maI

# Compute first day where BOTH:
#   D_total >= 2500
#   maI     >= 250
tau <- calc_trigger_tau_sird(sim$days, D_total, maI, CFG)

# Print initial fractions
print(sim$init_fraction)

# Print trigger summary
print_trigger_summary_sird(sim$days, D_total, maI, tau, CFG)

# Plot deaths and maI with trigger markers
plot_trigger_all_in_one_sird(sim$days, D_total, maI, tau, CFG)