############################################################
# Algorithm 1: maI(t) and maD(t) on ONE GRAPH (M = 1) — LOG y-axis
# FULL CORRECT VERSION
#
# What this does:
# 1) Simulate cumulative I(t), D(t) with correlated shocks (rho)
# 2) Convert to daily new: ΔI(t), ΔD(t)
# 3) Compute 7-day moving averages: maI(t), maD(t)
# 4) Plot maI and maD on ONE graph with log y-axis
#    - Uses your chosen y_min/y_max (must be > 0 for log scale)
#    - Adds threshold lines (theta_I, theta_D)
#    - Adds first/last crossing vertical lines for both series
############################################################

#############################
# 0) ONE-PLACE CONFIG
#############################
CFG <- list(
  # --- Simulation (Algorithm 1) ---
  N      = 1104,
  dt     = 1,
  M      = 1,          # M=1 for clarity
  I0     = 3,
  D0     = 1,
  KI     = 833e6,
  KD     = 0.4e6,
  gI     = 0.320,
  gD     = 0.320,
  sigmaI = 0.1,
  sigmaD = 0.1,
  rho    = 0.5,
  seed   = NULL,       # set 42 for reproducible; NULL for random each run
  
  # --- Thresholds (moving-average thresholds) ---
  theta_I = 5000,
  theta_D = 2500,
  
  # --- Plot window (x) ---
  x_min = 0,
  x_max = 150,         # set 1104 for full horizon
  
  # --- Plot window (y) ---
  y_min = 1,           # MUST be > 0 for log scale
  y_max = NA,          # NA => AUTO y_max; or set a number like 1e6
  
  # --- Log-safe floor for plotting-only ---
  eps_y = 1e-6
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
  rho <- cfg$rho; seed <- cfg$seed

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
# 2) CUMULATIVE -> DAILY NEW: ΔX(t) = max(0, X(t)-X(t-1))
############################################################
calc_daily_new <- function(X_mat) {
  dX <- X_mat
  dX[1, ] <- 0
  dX[-1, ] <- pmax(0, X_mat[-1, ] - X_mat[-nrow(X_mat), ])
  dX
}


############################################################
# 3) 7-DAY MOVING AVERAGE for a vector (M=1)
############################################################
calc_ma7_vec <- function(dX_vec) {
  as.numeric(stats::filter(dX_vec, rep(1/7, 7), sides = 1))
}


############################################################
# 4) FIRST & LAST crossing day for a single path x(t) >= threshold
############################################################
first_last_cross_single <- function(x, threshold_y) {
  idx <- which(x >= threshold_y)
  if (length(idx) == 0) return(list(first = NA_integer_, last = NA_integer_))
  list(first = idx[1] - 1, last = idx[length(idx)] - 1)
}


############################################################
# 5) PLOT maI and maD on ONE GRAPH with LOG y-axis
############################################################
plot_maI_maD_one_graph_logy <- function(days, maI, maD, cfg) {
  
  x_min <- cfg$x_min
  x_max <- cfg$x_max
  theta_I <- cfg$theta_I
  theta_D <- cfg$theta_D
  eps_y <- cfg$eps_y
  
  # ---- log-safe versions for plotting (replace NA/0 with eps) ----
  maI_plot <- maI
  maD_plot <- maD
  maI_plot[is.na(maI_plot) | maI_plot <= 0] <- eps_y
  maD_plot[is.na(maD_plot) | maD_plot <= 0] <- eps_y
  
  theta_I_plot <- max(theta_I, eps_y)
  theta_D_plot <- max(theta_D, eps_y)
  
  # ---- y limits: respect user y_min; auto y_max if NA ----
  y_min_plot <- max(cfg$y_min, eps_y)  # must be > 0 for log scale
  
  if (is.na(cfg$y_max)) {
    y_max_plot <- max(c(maI_plot, maD_plot, theta_I_plot, theta_D_plot), na.rm = TRUE) * 1.10
  } else {
    y_max_plot <- cfg$y_max
  }
  
  if (y_max_plot <= y_min_plot) y_max_plot <- y_min_plot * 10
  
  # ---- colors ----
  col_I <- "blue"
  col_D <- "orange"
  col_first <- "darkgreen"
  col_last  <- "purple"
  
  # ---- plot maI ----
  plot(
    x = days, y = maI_plot,
    type = "l", lwd = 2, col = col_I,
    log  = "y",
    xlab = "Days",
    ylab = "7-day moving average (daily new) [log scale]",
    main = "Algorithm 1: maI(t) and maD(t) (M = 1) — log y",
    xlim = c(x_min, x_max),
    ylim = c(y_min_plot, y_max_plot)
  )
  
  # ---- add maD ----
  lines(days, maD_plot, lwd = 2, col = col_D)
  
  # ---- thresholds ----
  abline(h = theta_I_plot, lwd = 2, lty = 3, col = col_I)
  abline(h = theta_D_plot, lwd = 2, lty = 3, col = col_D)
  
  text(x_min + 0.65 * (x_max - x_min), theta_I_plot,
       paste0("theta_I = ", format(theta_I, big.mark=",")),
       col = col_I, pos = 3)
  text(x_min + 0.65 * (x_max - x_min), theta_D_plot,
       paste0("theta_D = ", format(theta_D, big.mark=",")),
       col = col_D, pos = 3)
  
  # ---- crossings (use original maI/maD, not eps-adjusted) ----
  cross_I <- first_last_cross_single(maI, theta_I)
  cross_D <- first_last_cross_single(maD, theta_D)
  
  if (is.finite(cross_I$first)) {
    abline(v = cross_I$first, lwd = 2, lty = 2, col = col_first)
    text(cross_I$first, y_max_plot / 1.6, paste0("I first = ", cross_I$first), col = col_first, pos = 4)
  }
  if (is.finite(cross_I$last)) {
    abline(v = cross_I$last, lwd = 2, lty = 4, col = col_last)
    text(cross_I$last, y_max_plot / 2.2, paste0("I last = ", cross_I$last), col = col_last, pos = 4)
  }
  
  if (is.finite(cross_D$first)) {
    abline(v = cross_D$first, lwd = 2, lty = 2, col = col_first)
    text(cross_D$first, y_max_plot / 3.0, paste0("D first = ", cross_D$first), col = col_first, pos = 4)
  }
  if (is.finite(cross_D$last)) {
    abline(v = cross_D$last, lwd = 2, lty = 4, col = col_last)
    text(cross_D$last, y_max_plot / 4.2, paste0("D last = ", cross_D$last), col = col_last, pos = 4)
  }
  
  legend(
    "topleft",
    legend = c("maI(t)", "maD(t)", "theta_I", "theta_D", "first crossing", "last crossing"),
    col    = c(col_I, col_D, col_I, col_D, col_first, col_last),
    lty    = c(1, 1, 3, 3, 2, 4),
    lwd    = c(2, 2, 2, 2, 2, 2),
    bty    = "n"
  )
  
  invisible(list(cross_I = cross_I, cross_D = cross_D,
                 y_min_plot = y_min_plot, y_max_plot = y_max_plot))
}


#############################
# 6) RUN
#############################
sim <- simulate_alg1_ID(CFG)

dI <- as.numeric(calc_daily_new(sim$I)[, 1])
dD <- as.numeric(calc_daily_new(sim$D)[, 1])

maI <- calc_ma7_vec(dI)
maD <- calc_ma7_vec(dD)

out <- plot_maI_maD_one_graph_logy(sim$days, maI, maD, CFG)

cat("I crossings: first =", out$cross_I$first, ", last =", out$cross_I$last, "\n")
cat("D crossings: first =", out$cross_D$first, ", last =", out$cross_D$last, "\n")
cat("Plot y-limits used: [", out$y_min_plot, ", ", out$y_max_plot, "]\n")