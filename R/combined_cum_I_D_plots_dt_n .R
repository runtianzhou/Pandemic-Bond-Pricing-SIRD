############################################################
# CUMULATIVE I(t) & D(t) — Algorithm 1 settings (ONE GRAPH, M = 1)
# - simulate I and D with correlation rho
# - plot both curves on one diagram
# - different colors for I vs D, thresholds, and crossing lines
#
# ADDITION:
# - dt_n internal step size
# - print how many Euler–Maruyama steps actually executed
############################################################

#############################
# 0) ONE-PLACE CONFIG
#############################
CFG <- list(
  # Simulation parameters
  N      = 1104,
  dt     = 1,           # 1 day on x-axis (daily storage)
  dt_n   = 0.1,       # internal sub-step size within each day
  M      = 1,           # number of Monte Carlo paths
  I0     = 3,
  D0     = 1,
  KI     = 833e6,
  KD     = 0.4e6,
  gI     = 0.320,
  gD     = 0.320,
  sigmaI = 0.1,
  sigmaD = 0.1,
  rho    = 0.5,
  seed   = 123,        # NULL => random each run
  
  threshold_I = 5000,
  threshold_D = 2500,
  
  x_min = 0,
  x_max = 100,
  y_min = 0,
  y_max = NA,
  
  use_log_y = TRUE
)


############################################################
# 1) SIMULATE CUMULATIVE I(t), D(t)
#    Euler–Maruyama + correlation
#    - store daily (dt=1)
#    - simulate internally using dt_n
############################################################
simulate_alg1_ID <- function(cfg) {
  
  N  <- cfg$N
  dt <- cfg$dt
  dt_n <- cfg$dt_n
  M <- cfg$M
  
  I0 <- cfg$I0
  D0 <- cfg$D0
  
  KI <- cfg$KI
  KD <- cfg$KD
  gI <- cfg$gI
  gD <- cfg$gD
  sigmaI <- cfg$sigmaI
  sigmaD <- cfg$sigmaD
  rho <- cfg$rho
  
  seed <- cfg$seed
  
  if (N < 1) stop("N must be >= 1")
  if (dt <= 0) stop("dt must be > 0")
  if (is.null(dt_n) || dt_n <= 0) stop("dt_n must be > 0")
  if (dt_n > dt) stop("dt_n must be <= dt (dt is one day)")
  if (M < 1) stop("M must be >= 1")
  if (KI <= 0 || KD <= 0) stop("KI and KD must be > 0")
  if (rho < -1 || rho > 1) stop("rho must be in [-1, 1]")
  
  # Make each day exactly n_sub sub-steps
  n_sub <- as.integer(round(dt / dt_n))
  if (n_sub < 1) n_sub <- 1
  dt_eff <- dt / n_sub      # effective time step actually used by the algorithm
  
  set.seed(seed)
  
  I <- matrix(NA_real_, nrow = N + 1, ncol = M)
  D <- matrix(NA_real_, nrow = N + 1, ncol = M)
  I[1, ] <- I0
  D[1, ] <- D0
  
  sqrt_dt_eff <- sqrt(dt_eff)
  sqrt_1mrho2 <- sqrt(1 - rho^2)
  
  # counters (for printing)
  em_steps_performed <- 0L
  normals_generated  <- 0L
  
  for (t_day in 2:(N + 1)) {
    
    I_curr <- I[t_day - 1, ] # the current value of I during the internal sub-stepping
    D_curr <- D[t_day - 1, ] # the current value of D
                                # t_day = outer loop (days)
    for (k in 1:n_sub) {        # k = inner loop (sub-steps inside a day)
      
      # each sub-step uses TWO independent normal vectors of length M
      dW1 <- rnorm(M, mean = 0, sd = sqrt_dt_eff)
      dW2 <- rnorm(M, mean = 0, sd = sqrt_dt_eff)
      
      dBI <- dW1
      dBD <- rho * dW1 + sqrt_1mrho2 * dW2
      
      I_new <- I_curr + I_curr * (1 - I_curr / KI) * (gI * dt_eff + sigmaI * dBI)
      D_new <- D_curr + D_curr * (1 - D_curr / KD) * (gD * dt_eff + sigmaD * dBD)
      
      I_curr <- pmin(pmax(I_new, 0), KI)
      D_curr <- pmin(pmax(D_new, 0), KD)
      
      em_steps_performed <- em_steps_performed + M
      normals_generated  <- normals_generated  + 2L * M
    }
    
    I[t_day, ] <- I_curr
    D[t_day, ] <- D_curr
  }
  
  list(
    I = I, D = D, days = 0:N,
    n_sub = n_sub, dt_eff = dt_eff,
    em_steps_performed = em_steps_performed,
    normals_generated = normals_generated
  )
}


############################################################
# 2) FIRST CROSSING DAY: first day X(t) >= threshold
############################################################
first_cross_day_single <- function(x_vec, threshold_y) {
  idx <- which(x_vec >= threshold_y)
  if (length(idx) == 0) return(NA_integer_)
  idx[1] - 1
}


############################################################
# 3) PLOT I and D ON ONE GRAPH (different colors)
############################################################
plot_cumulative_I_and_D_one_graph <- function(sim, cfg) {
  
  days  <- sim$days
  I_vec <- as.numeric(sim$I[, 1])
  D_vec <- as.numeric(sim$D[, 1])
  
  x_min <- cfg$x_min
  x_max <- cfg$x_max
  use_log_y <- cfg$use_log_y
  
  col_I <- "blue"
  col_D <- "orange"
  
  if (use_log_y) {
    y_min <- max(cfg$y_min, 1e-6)
    if (is.na(cfg$y_max)) {
      y_max <- max(c(I_vec, D_vec), na.rm = TRUE) * 1.05
    } else {
      y_max <- cfg$y_max
    }
    I_plot <- pmax(I_vec, 1e-6)
    D_plot <- pmax(D_vec, 1e-6)
  } else {
    y_min <- cfg$y_min
    if (is.na(cfg$y_max)) {
      y_max <- max(c(I_vec, D_vec), na.rm = TRUE) * 1.05
    } else {
      y_max <- cfg$y_max
    }
    I_plot <- I_vec
    D_plot <- D_vec
  }
  
  plot(
    x = days, y = I_plot,
    type = "l", lwd = 2, col = col_I,
    log  = if (use_log_y) "y" else "",
    xlab = "Days",
    ylab = if (use_log_y) "Cumulative counts (log scale)" else "Cumulative counts",
    main = "Algorithm 1: Cumulative I(t) and D(t) (M = 1) − log y",
    xlim = c(x_min, x_max),
    ylim = c(y_min, y_max)
  )
  
  lines(days, D_plot, lwd = 2, col = col_D)
  
  legend(
    "topleft",
    legend = c("I(t)", "D(t)"),
    col    = c(col_I, col_D),
    lty    = c(1, 1),
    lwd    = c(2, 2),
    bty    = "n"
  )
}



# 
# plot_cumulative_I_and_D_one_graph <- function(sim, cfg) {
#   
#   days  <- sim$days
#   I_vec <- as.numeric(sim$I[, 1])
#   D_vec <- as.numeric(sim$D[, 1])
#   
#   x_min <- cfg$x_min
#   x_max <- cfg$x_max
#   use_log_y <- cfg$use_log_y
#   
#   col_I <- "blue"
#   col_D <- "orange"
#   col_thrI <- "blue"
#   col_thrD <- "orange"
#   col_crossI <- "darkgreen"
#   col_crossD <- "purple"
#   
#   if (use_log_y) {
#     y_min <- max(cfg$y_min, 1e-6)
#     if (is.na(cfg$y_max)) {
#       y_max <- max(c(I_vec, D_vec, cfg$threshold_I, cfg$threshold_D), na.rm = TRUE) * 1.05
#     } else {
#       y_max <- cfg$y_max
#     }
#     I_plot <- pmax(I_vec, 1e-6)
#     D_plot <- pmax(D_vec, 1e-6)
#   } else {
#     y_min <- cfg$y_min
#     if (is.na(cfg$y_max)) {
#       y_max <- max(c(I_vec, D_vec, cfg$threshold_I, cfg$threshold_D), na.rm = TRUE) * 1.05
#     } else {
#       y_max <- cfg$y_max
#     }
#     I_plot <- I_vec
#     D_plot <- D_vec
#   }
#   
#   plot(
#     x = days, y = I_plot,
#     type = "l", lwd = 2, col = col_I,
#     log  = if (use_log_y) "y" else "",
#     xlab = "Days",
#     ylab = if (use_log_y) "Cumulative counts (log scale)" else "Cumulative counts",
#     main = "Algorithm 1: Cumulative I(t) and D(t) (daily output, sub-stepped simulation)",
#     xlim = c(x_min, x_max),
#     ylim = c(y_min, y_max)
#   )
#   
#   lines(days, D_plot, lwd = 2, col = col_D)
#   
#   abline(h = cfg$threshold_I, lwd = 2, lty = 3, col = col_thrI)
#   abline(h = cfg$threshold_D, lwd = 2, lty = 3, col = col_thrD)
# 
#   text(x_min + 0.65 * (x_max - x_min), cfg$threshold_I,
#        paste0("threshold_I = ", format(cfg$threshold_I, big.mark=",")),
#        col = col_thrI, pos = 3)
#   text(x_min + 0.65 * (x_max - x_min), cfg$threshold_D,
#        paste0("threshold_D = ", format(cfg$threshold_D, big.mark=",")),
#        col = col_thrD, pos = 3)
# 
#   tau_I <- first_cross_day_single(I_vec, cfg$threshold_I)
#   tau_D <- first_cross_day_single(D_vec, cfg$threshold_D)
# 
#   if (is.finite(tau_I)) {
#     abline(v = tau_I, lwd = 2, lty = 2, col = col_crossI)
#     text(tau_I, y_max / 1.4, paste0("I crosses at day ", tau_I), col = col_crossI, pos = 4)
#   }
#   if (is.finite(tau_D)) {
#     abline(v = tau_D, lwd = 2, lty = 2, col = col_crossD)
#     text(tau_D, y_max / 2.0, paste0("D crosses at day ", tau_D), col = col_crossD, pos = 4)
#   }
#   
#   legend(
#     "topleft",
#     legend = c("I(t)", "D(t)", "threshold_I", "threshold_D", "I crossing day", "D crossing day"),
#     col    = c(col_I, col_D, col_thrI, col_thrD, col_crossI, col_crossD),
#     lty    = c(1, 1, 3, 3, 2, 2),
#     lwd    = c(2, 2, 2, 2, 2, 2),
#     bty    = "n"
#   )
#   
#   invisible(list(tau_I = tau_I, tau_D = tau_D, y_min = y_min, y_max = y_max))
# }


############################################################
# 4) RUN
############################################################
sim <- simulate_alg1_ID(CFG)
out <- plot_cumulative_I_and_D_one_graph(sim, CFG)

cat("========================================\n")
cat("SIMULATION SUMMARY\n")
cat("========================================\n")
cat("Paths (M)                    =", CFG$M, "\n")
cat("Days simulated (N)           =", CFG$N, "\n")
cat("Daily storage step (dt)      =", CFG$dt, "\n")
cat("Requested internal step dt_n =", CFG$dt_n, "\n")
cat("Internal sub-steps per day   =", sim$n_sub, "\n")
cat("Actual internal dt_eff       =", sim$dt_eff, "\n\n")

cat("Euler–Maruyama updates performed (all paths, all sub-steps) =",
    sim$em_steps_performed, "\n")
cat("Total rnorm() draws used (scalar normals) =",
    sim$normals_generated, "\n")
cat("========================================\n\n")

cat("I first crossing day =", out$tau_I, "\n")
cat("D first crossing day =", out$tau_D, "\n")