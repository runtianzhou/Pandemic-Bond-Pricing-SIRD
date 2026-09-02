############################################################
# CUMULATIVE I(t) & D(t) — Algorithm 1 settings (ONE GRAPH, M = 1)
# - simulate I and D with correlation rho
# - plot both curves on one diagram
# - different colors for I vs D, thresholds, and crossing lines
############################################################

#############################
# 0) ONE-PLACE CONFIG
#############################
CFG <- list(
  # Simulation parameters
  N      = 1104,
  dt     = 1,
  M      = 1,          # <-- M=1 for clear view
  I0     = 3,          # I_0 > D_0
  D0     = 1,
  KI     = 833e6,
  KD     = 0.4e6,
  gI     = 0.320,
  gD     = 0.320,
  sigmaI = 0.1,
  sigmaD = 0.1,
  rho    = 0.5,
  seed   = NULL,         # set NULL for random each run
  
  # Threshold lines (reference levels on cumulative plots)
  threshold_I = 5000,
  threshold_D = 2500,
  
  # Plot window
  x_min = 0,
  x_max = 100,         # set 1104 for full horizon
  y_min = 0,           # if using log="y", must be >0; if linear, you can set 0
  y_max = NA,          # NA => auto
  
  # Plot style
  use_log_y = TRUE      # TRUE for log-y to see both I and D; FALSE for linear
)


############################################################
# 1) SIMULATE CUMULATIVE I(t), D(t) (Euler–Maruyama + correlation)
############################################################
simulate_alg1_ID <- function(cfg) {
  N <- cfg$N; dt <- cfg$dt; M <- cfg$M
  I0 <- cfg$I0; D0 <- cfg$D0
  KI <- cfg$KI; KD <- cfg$KD
  gI <- cfg$gI; gD <- cfg$gD
  sigmaI <- cfg$sigmaI; sigmaD <- cfg$sigmaD
  rho <- cfg$rho
  seed <- cfg$seed

  if (N < 1) stop("N must be >= 1")
  if (dt <= 0) stop("dt must be > 0")
  if (M < 1) stop("M must be >= 1")
  if (KI <= 0 || KD <= 0) stop("KI and KD must be > 0")
  if (rho < -1 || rho > 1) stop("rho must be in [-1, 1]")

  set.seed(seed)

  I <- matrix(NA_real_, nrow = N + 1, ncol = M)
  D <- matrix(NA_real_, nrow = N + 1, ncol = M)
  I[1, ] <- I0
  D[1, ] <- D0

  sqrt_dt <- sqrt(dt)
  sqrt_1mrho2 <- sqrt(1 - rho^2)

  for (t in 2:(N + 1)) {
    dW1 <- rnorm(M, mean = 0, sd = sqrt_dt)
    dW2 <- rnorm(M, mean = 0, sd = sqrt_dt)

    dBI <- dW1
    dBD <- rho * dW1 + sqrt_1mrho2 * dW2

    I_prev <- I[t - 1, ]
    D_prev <- D[t - 1, ]

    I_new <- I_prev + I_prev * (1 - I_prev / KI) * (gI * dt + sigmaI * dBI)
    D_new <- D_prev + D_prev * (1 - D_prev / KD) * (gD * dt + sigmaD * dBD)

    I[t, ] <- pmin(pmax(I_new, 0), KI)
    D[t, ] <- pmin(pmax(D_new, 0), KD)
  }

  list(I = I, D = D, days = 0:N)
}

# 
# ############################################################
# # 1) SIMULATE CUMULATIVE I(t), D(t)
# #    Modified to match logistic2009.R model:
# #    dI = gI*I*(1-I/KI) dt + sigmaI*I*(1-I/KI) dBI
# #    dD = gD*D*(1-D/KD) dt + sigmaD*D*(1-D/KD) dBD
# #    with corr(dBI, dBD)=rho
# ############################################################
# simulate_alg1_ID <- function(cfg) {
#   N <- cfg$N; dt <- cfg$dt; M <- cfg$M
#   I0 <- cfg$I0; D0 <- cfg$D0
#   KI <- cfg$KI; KD <- cfg$KD
#   gI <- cfg$gI; gD <- cfg$gD
#   sigmaI <- cfg$sigmaI; sigmaD <- cfg$sigmaD
#   rho <- cfg$rho
#   seed <- cfg$seed
#   
#   if (N < 1) stop("N must be >= 1")
#   if (dt <= 0) stop("dt must be > 0")
#   if (M < 1) stop("M must be >= 1")
#   if (KI <= 0 || KD <= 0) stop("KI and KD must be > 0")
#   if (rho < -1 || rho > 1) stop("rho must be in [-1, 1]")
#   
#   set.seed(seed)
#   
#   I <- matrix(NA_real_, nrow = N + 1, ncol = M)
#   D <- matrix(NA_real_, nrow = N + 1, ncol = M)
#   I[1, ] <- I0
#   D[1, ] <- D0
#   
#   sqrt_dt <- sqrt(dt)
#   sqrt_1mrho2 <- sqrt(1 - rho^2)
#   
#   for (t in 2:(N + 1)) {
#     # Independent Brownian increments
#     dW1 <- rnorm(M, mean = 0, sd = sqrt_dt)
#     dW2 <- rnorm(M, mean = 0, sd = sqrt_dt)
#     
#     # Correlated increments matching logistic2009.R:
#     # dBI = dW1
#     # dBD = rho*dW1 + sqrt(1-rho^2)*dW2
#     dBI <- dW1
#     dBD <- rho * dW1 + sqrt_1mrho2 * dW2
#     
#     I_prev <- I[t - 1, ]
#     D_prev <- D[t - 1, ]
#     
#     # logistic2009.R drift and diffusion (Euler–Maruyama form)
#     drift_I <- gI     * I_prev * (1 - I_prev / KI)
#     diff_I  <- sigmaI * I_prev * (1 - I_prev / KI)
#     
#     drift_D <- gD     * D_prev * (1 - D_prev / KD)
#     diff_D  <- sigmaD * D_prev * (1 - D_prev / KD)
#     
#     I_new <- I_prev + drift_I * dt + diff_I * dBI
#     D_new <- D_prev + drift_D * dt + diff_D * dBD
#     
#     # Keep in admissible range (same idea as before)
#     I[t, ] <- pmin(pmax(I_new, 0), KI)
#     D[t, ] <- pmin(pmax(D_new, 0), KD)
#   }
#   
#   list(I = I, D = D, days = 0:N)
# }





############################################################
# 2) FIRST CROSSING DAY: first day X(t) >= threshold
############################################################
first_cross_day_single <- function(x_vec, threshold_y) {
  idx <- which(x_vec >= threshold_y)
  if (length(idx) == 0) return(NA_integer_)
  idx[1] - 1  # row1=day0
}


############################################################
# 3) PLOT I and D ON ONE GRAPH (different colors)
############################################################
plot_cumulative_I_and_D_one_graph <- function(sim, cfg) {
  
  days <- sim$days
  I_vec <- as.numeric(sim$I[, 1])
  D_vec <- as.numeric(sim$D[, 1])
  
  x_min <- cfg$x_min; x_max <- cfg$x_max
  use_log_y <- cfg$use_log_y
  
  # Colors
  col_I <- "blue"
  col_D <- "orange"
  col_thrI <- "blue"
  col_thrD <- "orange"
  col_crossI <- "darkgreen"
  col_crossD <- "purple"
  
  # y-limits
  if (use_log_y) {
    # log scale requires positive ylim
    y_min <- max(cfg$y_min, 1e-6)
    # auto y_max if NA
    if (is.na(cfg$y_max)) {
      y_max <- max(c(I_vec, D_vec, cfg$threshold_I, cfg$threshold_D), na.rm = TRUE) * 1.05
    } else {
      y_max <- cfg$y_max
    }
    # protect against zeros in data on log plot (for plotting only)
    I_plot <- pmax(I_vec, 1e-6)
    D_plot <- pmax(D_vec, 1e-6)
  } else {
    y_min <- cfg$y_min
    if (is.na(cfg$y_max)) {
      y_max <- max(c(I_vec, D_vec, cfg$threshold_I, cfg$threshold_D), na.rm = TRUE) * 1.05
    } else {
      y_max <- cfg$y_max
    }
    I_plot <- I_vec
    D_plot <- D_vec
  }
  
  # Base plot with I(t)
  plot(
    x = days, y = I_plot,
    type = "l", lwd = 2, col = col_I,
    log  = if (use_log_y) "y" else "",
    xlab = "Days",
    ylab = if (use_log_y) "Cumulative counts (log scale)" else "Cumulative counts",
    main = "Algorithm 1: Cumulative I(t) and D(t)",
    xlim = c(x_min, x_max),
    ylim = c(y_min, y_max)
  )
  
  # Add D(t)
  lines(days, D_plot, lwd = 2, col = col_D)
  
  # Threshold lines (same colors as curves)
  abline(h = cfg$threshold_I, lwd = 2, lty = 3, col = col_thrI)
  abline(h = cfg$threshold_D, lwd = 2, lty = 3, col = col_thrD)
  
  # Threshold labels
  text(x_min + 0.65 * (x_max - x_min), cfg$threshold_I,
       paste0("threshold_I = ", format(cfg$threshold_I, big.mark=",")),
       col = col_thrI, pos = 3)
  text(x_min + 0.65 * (x_max - x_min), cfg$threshold_D,
       paste0("threshold_D = ", format(cfg$threshold_D, big.mark=",")),
       col = col_thrD, pos = 3)
  
  # First crossing days (one path)
  tau_I <- first_cross_day_single(I_vec, cfg$threshold_I)
  tau_D <- first_cross_day_single(D_vec, cfg$threshold_D)
  
  if (is.finite(tau_I)) {
    abline(v = tau_I, lwd = 2, lty = 2, col = col_crossI)
    text(tau_I, y_max / 1.4, paste0("I crosses at day ", tau_I), col = col_crossI, pos = 4)
  }
  if (is.finite(tau_D)) {
    abline(v = tau_D, lwd = 2, lty = 2, col = col_crossD)
    text(tau_D, y_max / 2.0, paste0("D crosses at day ", tau_D), col = col_crossD, pos = 4)
  }
  
  legend(
    "topleft",
    legend = c("I(t)", "D(t)", "threshold_I", "threshold_D", "I crossing day", "D crossing day"),
    col    = c(col_I, col_D, col_thrI, col_thrD, col_crossI, col_crossD),
    lty    = c(1, 1, 3, 3, 2, 2),
    lwd    = c(2, 2, 2, 2, 2, 2),
    bty    = "n"
  )
  
  invisible(list(tau_I = tau_I, tau_D = tau_D, y_min = y_min, y_max = y_max))
}


############################################################
# 4) RUN
############################################################
sim <- simulate_alg1_ID(CFG)
out <- plot_cumulative_I_and_D_one_graph(sim, CFG)

cat("I first crossing day =", out$tau_I, "\n")
cat("D first crossing day =", out$tau_D, "\n")