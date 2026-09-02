############################################################
# ONE GRAPH: maI(t), maD(t), GR(t) (SE) all together
# - 3 y-axes (maI left, maD right, GR far-right)
# - all M paths shown
# - thresholds + trigger day lines
############################################################

#############################
# 0) CONFIG
#############################
CFG <- list(
  N  = 1104,
  dt = 1,
  M  = 1,        # set M
  I0 = 1,
  D0 = 1,
  KI = 833e6,
  KD = 0.4e6,
  gI = 0.320,
  gD = 0.320,
  sigmaI = 0.1,
  sigmaD = 0.1,
  rho    = 0.5,
  seed   = 123,
  
  theta_I = 5000,
  theta_D = 2500,
  T_maturity = 1104,
  
  x_min = 0,
  x_max = 200
)

#############################
# 1) SIMULATE I(t), D(t)
#############################
simulate_alg1_ID <- function(cfg) {
  N <- cfg$N; dt <- cfg$dt; M <- cfg$M
  I0 <- cfg$I0; D0 <- cfg$D0
  KI <- cfg$KI; KD <- cfg$KD
  gI <- cfg$gI; gD <- cfg$gD
  sigmaI <- cfg$sigmaI; sigmaD <- cfg$sigmaD
  rho <- cfg$rho; seed <- cfg$seed
  
  set.seed(seed)
  
  I <- matrix(NA_real_, nrow = N + 1, ncol = M)
  D <- matrix(NA_real_, nrow = N + 1, ncol = M)
  I[1, ] <- I0
  D[1, ] <- D0
  
  sqrt_dt <- sqrt(dt)
  sqrt_1mrho2 <- sqrt(1 - rho^2)
  
  for (t in 2:(N + 1)) {
    dW1 <- rnorm(M, 0, sqrt_dt)
    dW2 <- rnorm(M, 0, sqrt_dt)
    
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

#############################
# 2) DAILY NEW
#############################
calc_daily_new <- function(X_mat) {
  dX <- X_mat
  dX[1, ] <- 0
  dX[-1, ] <- pmax(0, X_mat[-1, ] - X_mat[-nrow(X_mat), ])
  dX
}

#############################
# 3) 7-day MA and SE
#############################
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

calc_GR_SE <- function(dI_mat) {
  maI <- ma7_mat(dI_mat)
  seI <- sd7_mat(dI_mat) / sqrt(7)     # SE version
  GR  <- maI - 1.533 * seI
  list(maI = maI, seI = seI, GR = GR)
}

#############################
# 4) Trigger tau per path
#############################
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
  tau
}

#############################
# 5) ONE GRAPH WITH 3 Y-AXES
#############################
plot_all_in_one <- function(days, maI, maD, GR, tau, cfg) {
  
  x_min <- cfg$x_min; x_max <- cfg$x_max
  M <- ncol(maI)
  
  # Colors (clear families)
  col_I <- "blue"
  col_D <- "orange"
  col_G <- "purple"
  
  # Build y-limits for each variable (avoid NA)
  yI_rng <- range(maI, na.rm = TRUE)
  yD_rng <- range(maD, na.rm = TRUE)
  yG_rng <- range(GR,  na.rm = TRUE)
  
  # --- Panel setup: extra right margin for 3rd axis ---
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mar = c(5, 5, 4, 6))  # extra right margin
  
  # ---- 1) Base plot: maI on LEFT axis ----
  matplot(days, maI, type="l", lty=1, lwd=2,
          col = col_I,
          xlab="Days", ylab="maI(t) (left axis)",
          main=paste0("Algorithm 1: maI (blue), maD (orange), GR (purple), M=", M),
          xlim=c(x_min, x_max),
          ylim=yI_rng)
  
  # maI threshold
  abline(h = cfg$theta_I, col = col_I, lty = 3, lwd = 2)
  
  # trigger lines (each path)
  tau_ok <- tau[is.finite(tau)]
  if (length(tau_ok) > 0) abline(v = tau_ok, col="red", lty=2, lwd=1)
  
  # ---- 2) Add maD on RIGHT axis ----
  par(new = TRUE)
  matplot(days, maD, type="l", lty=1, lwd=2,
          col = col_D,
          axes = FALSE, xlab="", ylab="",
          xlim=c(x_min, x_max),
          ylim=yD_rng)
  axis(side = 4, col = col_D, col.axis = col_D)
  mtext("maD(t) (right axis)", side = 4, line = 2, col = col_D)
  
  # maD threshold
  abline(h = cfg$theta_D, col = col_D, lty = 3, lwd = 2)
  
  # ---- 3) Add GR on FAR-RIGHT axis (shifted outward) ----
  par(new = TRUE)
  matplot(days, GR, type="l", lty=1, lwd=2,
          col = col_G,
          axes = FALSE, xlab="", ylab="",
          xlim=c(x_min, x_max),
          ylim=yG_rng)
  axis(side = 4, line = 3.7, col = col_G, col.axis = col_G)
  mtext("GR(t) (far-right axis)", side = 4, line = 5.5, col = col_G)
  
  # GR=0 line
  abline(h = 0, col = col_G, lty = 2, lwd = 2)
  
  # Legend
  legend("topleft",
         legend = c("maI paths", "maD paths", "GR paths", "theta_I", "theta_D", "GR=0", "trigger τ (each path)"),
         col    = c(col_I, col_D, col_G, col_I, col_D, col_G, "red"),
         lty    = c(1, 1, 1, 3, 3, 2, 2),
         lwd    = c(2, 2, 2, 2, 2, 2, 1),
         bty    = "n")
}

#############################
# RUN
#############################
sim <- simulate_alg1_ID(CFG)

dI <- calc_daily_new(sim$I)
dD <- calc_daily_new(sim$D)

gr_out <- calc_GR_SE(dI)
maI <- gr_out$maI
GR  <- gr_out$GR
maD <- ma7_mat(dD)

tau <- calc_trigger_tau(maI, maD, GR, CFG)
cat("tau per path:\n"); print(tau)

plot_all_in_one(sim$days, maI, maD, GR, tau, CFG)