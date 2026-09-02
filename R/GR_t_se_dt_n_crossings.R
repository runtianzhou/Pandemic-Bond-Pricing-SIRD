############################################################
# GR(t) (SE version) + intervals GR(t)>0
# FULL CORRECT VERSION (WITH dt_n INTERNAL SUB-STEPPING)
#
# Keep dt = 1 as "one day" for x-axis + storage,
# add dt_n to refine Euler–Maruyama simulation inside each day.
############################################################

#############################
# 0) CONFIG
#############################
CFG <- list(
  # Simulation
  N      = 1104,
  dt     = 1,            # <- 1 day on x-axis (daily storage)
  dt_n   = 0.1,        # <- NEW: internal sub-step within each day (try 0.0001)
  M      = 1,            # small (5–20) for readable plot/labels
  I0     = 1,
  D0     = 1,
  KI     = 833e6,
  KD     = 0.4e6,
  gI     = 0.320,
  gD     = 0.320,
  sigmaI = 0.1,
  sigmaD = 0.1,
  rho    = 0.5,
  seed   = 123,
  
  # Plot window
  x_min = 0,
  x_max = 200,
  y_min = NA,
  y_max = NA,
  
  # Style
  lwd_gr  = 2,
  lwd_v   = 2,
  cex_lab = 0.75,
  
  # Output
  write_csv = FALSE,
  csv_name  = "GR_positive_intervals.csv"
)

############################################################
# 1) SIMULATE I(t), D(t) with dt_n sub-stepping
############################################################
simulate_alg1_ID <- function(cfg) {
  
  N  <- cfg$N
  dt <- cfg$dt
  dt_n <- cfg$dt_n
  M  <- cfg$M
  
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
  
  # Counters (optional but nice)
  em_steps_performed <- 0L
  normals_generated  <- 0L
  
  for (t_day in 2:(N + 1)) {
    
    I_curr <- I[t_day - 1, ]
    D_curr <- D[t_day - 1, ]
    
    for (k in 1:n_sub) {
      
      dW1 <- rnorm(M, 0, sqrt_dt_eff)
      dW2 <- rnorm(M, 0, sqrt_dt_eff)
      
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
# 2) Daily new
############################################################
calc_daily_new <- function(X_mat) {
  dX <- X_mat
  dX[1, ] <- 0
  dX[-1, ] <- pmax(0, X_mat[-1, ] - X_mat[-nrow(X_mat), ])
  dX
}

############################################################
# 3) 7-day MA, SD, SE
############################################################
ma7_mat <- function(dX_mat) {
  n <- nrow(dX_mat); M <- ncol(dX_mat)
  ma <- matrix(NA_real_, nrow = n, ncol = M)
  for (j in 1:M) ma[, j] <- stats::filter(dX_mat[, j], rep(1/7, 7), sides = 1)
  ma
}

sd7_mat <- function(dX_mat) {
  n <- nrow(dX_mat); M <- ncol(dX_mat)
  out <- matrix(NA_real_, nrow = n, ncol = M)
  for (j in 1:M) {
    x <- dX_mat[, j]
    for (t in 7:n) out[t, j] <- sd(x[(t-6):t])
  }
  out
}

calc_GR_mat_SE <- function(dI_mat) {
  maI <- ma7_mat(dI_mat)
  sdI <- sd7_mat(dI_mat)
  seI <- sdI / sqrt(7)              # SE version
  GR  <- maI - 1.533 * seI
  list(maI = maI, seI = seI, GR = GR)
}

############################################################
# 4) GR>0 intervals for ONE path (vector)
############################################################
gr_positive_intervals_vec <- function(GR_vec) {
  pos <- (!is.na(GR_vec)) & (GR_vec > 0)
  
  if (!any(pos)) {
    return(data.frame(start_day=integer(0), end_day=integer(0), duration=integer(0)))
  }
  
  d <- diff(c(FALSE, pos, FALSE))   # +1 at start, -1 after end
  starts <- which(d == 1)
  ends   <- which(d == -1) - 1
  
  start_day <- starts - 1           # row 1 = day 0
  end_day   <- ends - 1
  duration  <- end_day - start_day
  
  data.frame(start_day=start_day, end_day=end_day, duration=duration)
}

############################################################
# 5) GR>0 intervals for ALL paths (matrix) -> one big table
############################################################
gr_positive_intervals_mat <- function(GR_mat) {
  M <- ncol(GR_mat)
  out <- vector("list", M)
  keep <- 0L
  
  for (j in 1:M) {
    df <- gr_positive_intervals_vec(GR_mat[, j])
    if (nrow(df) > 0) {
      df$path <- j
      keep <- keep + 1L
      out[[keep]] <- df
    }
  }
  
  if (keep == 0L) {
    return(data.frame(path=integer(0), start_day=integer(0), end_day=integer(0), duration=integer(0)))
  }
  
  res <- do.call(rbind, out[1:keep])
  res <- res[, c("path", "start_day", "end_day", "duration")]
  rownames(res) <- NULL
  res
}

############################################################
# 6) Per-path summary: total duration, max duration, #segments
############################################################
summarize_intervals_by_path <- function(intervals_df, M) {
  total_dur <- integer(M)
  max_dur   <- integer(M)
  n_seg     <- integer(M)
  
  if (nrow(intervals_df) > 0) {
    for (j in 1:M) {
      dj <- intervals_df[intervals_df$path == j, , drop=FALSE]
      if (nrow(dj) > 0) {
        total_dur[j] <- sum(dj$duration)
        max_dur[j]   <- max(dj$duration)
        n_seg[j]     <- nrow(dj)
      }
    }
  }
  
  data.frame(path=1:M, segments=n_seg, total_duration=total_dur, max_segment=max_dur)
}

############################################################
# 7) Plot all GR paths + start/end lines + labels
############################################################
plot_GR_with_intervals <- function(days, GR_mat, intervals_df, cfg) {
  x_min <- cfg$x_min; x_max <- cfg$x_max
  y_min <- cfg$y_min; y_max <- cfg$y_max
  M <- ncol(GR_mat)
  
  if (is.na(y_min)) y_min <- min(GR_mat, na.rm = TRUE)
  if (is.na(y_max)) y_max <- max(GR_mat, na.rm = TRUE) * 1.05
  
  path_cols <- 1:M
  
  matplot(
    x = days, y = GR_mat,
    type="l", lty=1, lwd=cfg$lwd_gr, col=path_cols,
    xlab="Days", ylab="GR(t)",
    main=paste0("GR(t) (SE version), M=", M, " with GR>0 intervals"),
    xlim=c(x_min, x_max),
    ylim=c(y_min, y_max)
  )
  abline(h=0, lty=2)
  
  # y-levels for labels (stagger)
  y_levels <- seq(y_max*0.95, y_max*0.65, length.out = M)
  
  # draw start/end lines and labels for each interval
  for (j in 1:M) {
    dj <- intervals_df[intervals_df$path == j, , drop=FALSE]
    if (nrow(dj) == 0) next
    
    for (k in 1:nrow(dj)) {
      s <- dj$start_day[k]
      e <- dj$end_day[k]
      dur <- dj$duration[k]
      
      abline(v = s, col = path_cols[j], lwd = cfg$lwd_v, lty = 2)
      abline(v = e, col = path_cols[j], lwd = cfg$lwd_v, lty = 2)
      
      text(s, y_levels[j],
           paste0("j", j, " start=", s),
           col = path_cols[j], pos = 4, cex = cfg$cex_lab)
      text(e, y_levels[j] - 0.03*(y_max-y_min),
           paste0("j", j, " end=", e, " (dur=", dur, ")"),
           col = path_cols[j], pos = 4, cex = cfg$cex_lab)
    }
  }
  
  legend("top", legend=paste0("Path ", 1:M), col=path_cols, lty=1, lwd=cfg$lwd_gr, bty="n")
}

#############################
# RUN
#############################
sim <- simulate_alg1_ID(CFG)
dI  <- calc_daily_new(sim$I)
gr  <- calc_GR_mat_SE(dI)

intervals_all <- gr_positive_intervals_mat(gr$GR)

cat("\n========================================\n")
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

cat("All GR(t)>0 intervals (all paths):\n")
print(intervals_all)

summary_by_path <- summarize_intervals_by_path(intervals_all, CFG$M)
cat("\nSummary by path:\n")
print(summary_by_path)

plot_GR_with_intervals(sim$days, gr$GR, intervals_all, CFG)

if (isTRUE(CFG$write_csv)) {
  write.csv(intervals_all, CFG$csv_name, row.names = FALSE)
  cat("\nWrote interval table to:", CFG$csv_name, "\n")
}