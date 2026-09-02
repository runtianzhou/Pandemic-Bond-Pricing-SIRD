############################################################
# Algorithm 1 Trigger (SE version) + ALL PATHS plots
# + thresholds on maI/maD
# + trigger tau
# + GR panel: begin/end of main GR>0 mountain per path
# + user-controlled y-limits (min/max) for each panel
#
# FULL CORRECT VERSION (WITH dt_n INTERNAL SUB-STEPPING)
#
# Keep dt = 1 day for x-axis + storage,
# add dt_n (sub-step) to refine Euler–Maruyama simulation inside each day.
############################################################

#############################
# 0) ONE-PLACE CONFIG
#############################
CFG <- list(
  # Simulation horizon
  N  = 1104,
  dt = 1,            # <- daily storage and x-axis
  dt_n = 1,      # <- NEW: internal sub-step within each day (try 0.0001)
  M  = 5,            # set 1..20
  
  # Initial values
  I0 = 3,
  D0 = 1,
  
  # Logistic parameters
  KI = 833e6,
  KD = 0.4e6,
  gI = 0.320,
  gD = 0.320,
  sigmaI = 0.1,
  sigmaD = 0.1,
  rho    = 0.5,
  seed   = 123,
  
  # Trigger thresholds
  theta_I = 5000,
  theta_D = 2500,
  T_maturity = 1104,
  
  # Plot window (x)
  x_min = 5,
  x_max = 165,
  
  # Plot window (y) — set to NA for auto
  y_maI_min = 1000,
  y_maI_max = NA,
  y_maD_min = 1,
  y_maD_max = NA,
  y_GR_min  = NA,
  y_GR_max  = NA,
  
  # Labels/toggles
  label_each_tau = TRUE,
  label_thresholds = TRUE,
  label_GR_mountain = TRUE
)

############################################################
# 1) SIMULATE cumulative I(t), D(t) with dt_n sub-stepping
############################################################
simulate_alg1_ID <- function(cfg) {
  N <- cfg$N; dt <- cfg$dt; dt_n <- cfg$dt_n; M <- cfg$M
  I0 <- cfg$I0; D0 <- cfg$D0
  KI <- cfg$KI; KD <- cfg$KD
  gI <- cfg$gI; gD <- cfg$gD
  sigmaI <- cfg$sigmaI; sigmaD <- cfg$sigmaD
  rho <- cfg$rho; seed <- cfg$seed
  
  if (N < 1) stop("N must be >= 1")
  if (dt <= 0) stop("dt must be > 0")
  if (is.null(dt_n) || dt_n <= 0) stop("dt_n must be > 0")
  if (dt_n > dt) stop("dt_n must be <= dt (dt is one day)")
  if (M < 1) stop("M must be >= 1")
  if (KI <= 0 || KD <= 0) stop("KI and KD must be > 0")
  if (rho < -1 || rho > 1) stop("rho must be in [-1, 1]")
  
  # Partition each day into integer sub-steps so total = 1 day exactly
  n_sub  <- as.integer(round(dt / dt_n))
  if (n_sub < 1) n_sub <- 1
  dt_eff <- dt / n_sub
  
  set.seed(seed)
  
  I <- matrix(NA_real_, nrow = N + 1, ncol = M)
  D <- matrix(NA_real_, nrow = N + 1, ncol = M)
  I[1, ] <- I0
  D[1, ] <- D0
  
  sqrt_dt_eff <- sqrt(dt_eff)
  sqrt_1mrho2 <- sqrt(1 - rho^2)
  
  # Counters (optional)
  em_steps_performed <- 0L
  normals_generated  <- 0L
  
  # Day loop (store only end-of-day)
  for (t_day in 2:(N + 1)) {
    
    I_curr <- I[t_day - 1, ]
    D_curr <- D[t_day - 1, ]
    
    # Sub-step loop (dt_eff inside each day)
    for (k in 1:n_sub) {
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
    normals_generated  = normals_generated
  )
}

############################################################
# 2) DAILY NEW: ΔX(t)=max(0,X(t)-X(t-1))
############################################################
calc_daily_new <- function(X_mat) {
  dX <- X_mat
  dX[1, ] <- 0
  dX[-1, ] <- pmax(0, X_mat[-1, ] - X_mat[-nrow(X_mat), ])
  dX
}

############################################################
# 3) 7-day moving average
############################################################
ma7_mat <- function(dX_mat) {
  n <- nrow(dX_mat); M <- ncol(dX_mat)
  ma <- matrix(NA_real_, nrow = n, ncol = M)
  for (j in 1:M) ma[, j] <- stats::filter(dX_mat[, j], rep(1/7, 7), sides = 1)
  ma
}

############################################################
# 4) rolling SD and SE (SE = SD/sqrt(7))
############################################################
sd7_mat <- function(dX_mat) {
  n <- nrow(dX_mat); M <- ncol(dX_mat)
  out <- matrix(NA_real_, nrow = n, ncol = M)
  for (j in 1:M) {
    x <- dX_mat[, j]
    for (t in 7:n) out[t, j] <- sd(x[(t-6):t])
  }
  out
}
se7_mat <- function(dX_mat) sd7_mat(dX_mat) / sqrt(7)

############################################################
# 5) GR(t) = maI - 1.533*seI  (SE version)
############################################################
calc_GR_SE <- function(dI_mat) {
  maI <- ma7_mat(dI_mat)
  seI <- se7_mat(dI_mat)
  GR  <- maI - 1.533 * seI
  list(maI = maI, seI = seI, GR = GR)
}

############################################################
# 6) Trigger day tau per path
############################################################
calc_trigger_tau <- function(maI, maD, GR, cfg) {
  days <- 0:(nrow(maI) - 1)
  Tm <- min(cfg$T_maturity, max(days))
  M <- ncol(maI)
  
  tau <- rep(NA_integer_, M)
  for (j in 1:M) {
    cond <- (days <= Tm) &
      (!is.na(maI[, j])) & (!is.na(maD[, j])) & (!is.na(GR[, j])) &
      (maD[, j] > cfg$theta_D) &
      (maI[, j] > cfg$theta_I) &
      (GR[, j]  > 0)
    hit <- which(cond)
    if (length(hit) > 0) tau[j] <- days[hit[1]]
  }
  
  list(
    tau = tau,
    min_tau = if (all(is.na(tau))) NA_integer_ else min(tau, na.rm=TRUE),
    max_tau = if (all(is.na(tau))) NA_integer_ else max(tau, na.rm=TRUE)
  )
}

############################################################
# 7) MAIN GR>0 "mountain" interval per path (longest run where GR>0)
############################################################
gr_main_positive_interval_vec <- function(GR_vec) {
  pos <- (!is.na(GR_vec)) & (GR_vec > 0)
  r <- rle(pos)
  if (!any(r$values)) return(list(start=NA_integer_, end=NA_integer_, duration=0L))
  
  true_runs <- which(r$values)
  best_run  <- true_runs[which.max(r$lengths[true_runs])]
  
  end_pos   <- cumsum(r$lengths)[best_run]
  start_pos <- end_pos - r$lengths[best_run] + 1
  
  list(
    start = start_pos - 1,
    end   = end_pos - 1,
    duration = (end_pos - 1) - (start_pos - 1)
  )
}

############################################################
# 8) Plot ALL paths in each panel + thresholds + trigger + GR mountain markers
############################################################
plot_trigger_panels_all_paths <- function(days, maI, maD, GR, trig, cfg) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mfrow = c(3, 1), mar = c(4, 5, 3, 2))
  
  x_min <- cfg$x_min
  x_max <- cfg$x_max
  M <- ncol(maI)
  path_cols <- 1:M
  x_lab <- x_min + 0.65 * (x_max - x_min)
  
  # =========================================================
  # Panel 1: maI (log scale)
  # =========================================================
  maI_plot <- maI
  maI_plot[!is.finite(maI_plot) | maI_plot <= 0] <- NA
  
  if (all(is.na(maI_plot))) {
    y1_min <- 1
    y1_max <- 10
  } else {
    y1_min <- if (is.na(cfg$y_maI_min)) min(maI_plot, na.rm = TRUE) else cfg$y_maI_min
    y1_max <- if (is.na(cfg$y_maI_max)) max(maI_plot, na.rm = TRUE) * 1.05 else cfg$y_maI_max
    
    if (!is.finite(y1_min) || y1_min <= 0) {
      y1_min <- min(maI_plot[maI_plot > 0], na.rm = TRUE)
    }
    if (!is.finite(y1_max) || y1_max <= y1_min) {
      y1_max <- y1_min * 10
    }
  }
  
  matplot(days, maI_plot, type = "l", lty = 1, lwd = 2, col = path_cols,
          xlab = "Days", ylab = "maI(t)",
          main = paste0("maI(t) (7-day avg new infections), M=", M),
          xlim = c(x_min, x_max), ylim = c(y1_min, y1_max), log = "y")
  
  if (is.finite(cfg$theta_I) && cfg$theta_I > 0) {
    abline(h = cfg$theta_I, lty = 3, lwd = 2, col = "darkred")
    if (isTRUE(cfg$label_thresholds)) {
      text(x_lab, cfg$theta_I, paste0("theta_I=", cfg$theta_I),
           col = "darkred", pos = 3, cex = 1)
    }
  }
  
  if (is.finite(trig$min_tau)) abline(v = trig$min_tau, col = "red", lwd = 2, lty = 2)
  if (is.finite(trig$max_tau)) abline(v = trig$max_tau, col = "red", lwd = 2, lty = 3)
  
  
  # =========================================================
  # Panel 2: maD (log scale)
  # =========================================================
  maD_plot <- maD
  maD_plot[!is.finite(maD_plot) | maD_plot <= 0] <- NA
  
  if (all(is.na(maD_plot))) {
    y2_min <- 1
    y2_max <- 10
  } else {
    y2_min <- if (is.na(cfg$y_maD_min)) min(maD_plot, na.rm = TRUE) else cfg$y_maD_min
    y2_max <- if (is.na(cfg$y_maD_max)) max(maD_plot, na.rm = TRUE) * 1.05 else cfg$y_maD_max
    
    if (!is.finite(y2_min) || y2_min <= 0) {
      y2_min <- min(maD_plot[maD_plot > 0], na.rm = TRUE)
    }
    if (!is.finite(y2_max) || y2_max <= y2_min) {
      y2_max <- y2_min * 10
    }
  }
  
  matplot(days, maD_plot, type = "l", lty = 1, lwd = 2, col = path_cols,
          xlab = "Days", ylab = "maD(t)",
          main = paste0("maD(t) (7-day avg new deaths), M=", M),
          xlim = c(x_min, x_max), ylim = c(y2_min, y2_max), log = "y")
  
  if (is.finite(cfg$theta_D) && cfg$theta_D > 0) {
    abline(h = cfg$theta_D, lty = 3, lwd = 2, col = "darkred")
    if (isTRUE(cfg$label_thresholds)) {
      text(x_lab, cfg$theta_D, paste0("theta_D=", cfg$theta_D),
           col = "darkred", pos = 3, cex = 1)
    }
  }
  
  if (is.finite(trig$min_tau)) abline(v = trig$min_tau, col = "red", lwd = 2, lty = 2)
  if (is.finite(trig$max_tau)) abline(v = trig$max_tau, col = "red", lwd = 2, lty = 3)
  
  
  # =========================================================
  # Panel 3: GR (log scale)
  # =========================================================
  GR_plot <- GR
  GR_plot[!is.finite(GR_plot) | GR_plot <= 0] <- NA
  
  if (all(is.na(GR_plot))) {
    y3_min <- 1
    y3_max <- 10
  } else {
    y3_min <- if (is.na(cfg$y_GR_min)) min(GR_plot, na.rm = TRUE) else cfg$y_GR_min
    y3_max <- if (is.na(cfg$y_GR_max)) max(GR_plot, na.rm = TRUE) * 1.05 else cfg$y_GR_max
    
    if (!is.finite(y3_min) || y3_min <= 0) {
      y3_min <- min(GR_plot[GR_plot > 0], na.rm = TRUE)
    }
    if (!is.finite(y3_max) || y3_max <= y3_min) {
      y3_max <- y3_min * 10
    }
  }
  
  matplot(days, GR_plot, type = "l", lty = 1, lwd = 2, col = path_cols,
          xlab = "Days", ylab = "GR(t)",
          main = paste0("GR(t) = maI - 1.533*seI (SE), M=", M),
          xlim = c(x_min, x_max), ylim = c(y3_min, y3_max), log = "y")
  
  if (is.finite(trig$min_tau)) abline(v = trig$min_tau, col = "red", lwd = 2, lty = 2)
  if (is.finite(trig$max_tau)) abline(v = trig$max_tau, col = "red", lwd = 2, lty = 3)
  
  # For each path: mark main GR>0 mountain start/end
  y_levels <- exp(seq(log(y3_max * 0.95), log(y3_min * 1e12), length.out = M))
  for (j in 1:M) {
    interval <- gr_main_positive_interval_vec(GR[, j])
    if (!is.finite(interval$start)) next
    
    abline(v = interval$start, col = path_cols[j], lwd = 2, lty = 2)
    abline(v = interval$end,   col = path_cols[j], lwd = 2, lty = 2)
    
    if (isTRUE(cfg$label_GR_mountain)) {
      text(interval$start, y_levels[j],
           paste0("j", j, " start=", interval$start),
           col = path_cols[j], pos = 4, cex = 1)
      text(interval$end, y_levels[j] / 1.2,
           paste0("j", j, " end=", interval$end),
           col = path_cols[j], pos = 2, cex = 1)
    }
  }
  
  # Optional: label each tau_j
  if (isTRUE(cfg$label_each_tau)) {
    tau_y <- exp(seq(log(y3_max), log(y3_min), length.out = M + 2))[2:(M + 1)]
    for (j in 1:M) {
      if (is.finite(trig$tau[j])) {
        text(trig$tau[j], tau_y[j],
             paste0("j", j, " trigger_date=", trig$tau[j]),
             col = path_cols[j], pos = 4, cex = 1) #path_cols[j]
      }
    }
  }
  
  legend("bottomleft", legend = paste0("Path ", 1:M),
         col = path_cols, lty = 1, lwd = 2, bty = "n")
}

#############################
# RUN
#############################
sim <- simulate_alg1_ID(CFG)
dI  <- calc_daily_new(sim$I)
dD  <- calc_daily_new(sim$D)

gr_out <- calc_GR_SE(dI)
maI <- gr_out$maI
GR  <- gr_out$GR
maD <- ma7_mat(dD)

trig <- calc_trigger_tau(maI, maD, GR, CFG)

cat("========================================\n")
cat("SIMULATION SUMMARY\n")
cat("========================================\n")
cat("Paths (M)                    =", CFG$M, "\n")
cat("Days simulated (N)           =", CFG$N, "\n")
cat("Daily storage step (dt)      =", CFG$dt, "\n")
cat("Requested internal step dt_n =", CFG$dt_n, "\n")
cat("Internal sub-steps per day   =", sim$n_sub, "\n")
cat("Actual internal dt_eff       =", sim$dt_eff, "\n")
cat("Euler–Maruyama updates performed (total across all paths) =",
    sim$em_steps_performed, "\n")
cat("Total rnorm() draws used (scalar normals) =",
    sim$normals_generated, "\n")
cat("========================================\n\n")

cat("tau per path:\n")
print(trig$tau)

plot_trigger_panels_all_paths(sim$days, maI, maD, GR, trig, CFG)