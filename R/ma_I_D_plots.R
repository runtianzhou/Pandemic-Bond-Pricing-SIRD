############################################################
# Algorithm 1: maI(t) and maD(t) in ONE FILE
# - simulate cumulative I(t), D(t) with correlation rho
# - daily new ΔI, ΔD
# - 7-day moving averages maI, maD
# - plot maI + maD (side-by-side) with thresholds + first/last crossings
############################################################

#############################
# 0) ONE-PLACE CONFIG
#############################
CFG <- list(
  # ---- Simulation (Algorithm 1) ----
  N      = 1104,
  dt     = 1,
  M      = 1,        # change once here
  I0     = 3,
  D0     = 1,
  KI     = 833e6,
  KD     = 0.4e6,
  gI     = 0.320,
  gD     = 0.320,
  sigmaI = 0.1,
  sigmaD = 0.1,
  rho    = 0.5,
  seed   = 123,       # set 42 for reproducible
  
  # ---- Thresholds for moving averages ----
  theta_I = 5000,
  theta_D = 2500,
  
  # ---- Plot zoom controls (maI) ----
  x_min_I = 30,
  x_max_I = 110,
  y_min_I = 0,
  y_max_I = NA,        # NA => auto
  
  # ---- Plot zoom controls (maD) ----
  x_min_D = 30,
  x_max_D = 65,
  y_min_D = 0,
  y_max_D = NA,        # NA => auto
  
  # ---- Plot performance ----
  max_paths_to_draw = 5000
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
    
    # Two independent Brownian increments
    dW1 <- rnorm(M, mean = 0, sd = sqrt_dt)
    dW2 <- rnorm(M, mean = 0, sd = sqrt_dt)
    
    # Correlated increments (Algorithm 1)
    dBI <- dW1
    dBD <- rho * dW1 + sqrt_1mrho2 * dW2
    
    I_prev <- I[t - 1, ]
    D_prev <- D[t - 1, ]
    
    # Euler–Maruyama stochastic logistic updates
    I_new <- I_prev + I_prev * (1 - I_prev / KI) * (gI * dt + sigmaI * dBI)
    D_new <- D_prev + D_prev * (1 - D_prev / KD) * (gD * dt + sigmaD * dBD)
    
    # Clamp to [0, K]
    I[t, ] <- pmin(pmax(I_new, 0), KI)
    D[t, ] <- pmin(pmax(D_new, 0), KD)
  }
  
  list(I = I, D = D, days = 0:N)
}


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
# 3) 7-DAY MOVING AVERAGE (trailing): ma(t)=mean(Δ(t-6..t))
############################################################
calc_ma7 <- function(dX_mat) {
  n <- nrow(dX_mat)
  M <- ncol(dX_mat)
  ma <- matrix(NA_real_, nrow = n, ncol = M)
  for (j in 1:M) {
    ma[, j] <- stats::filter(dX_mat[, j], rep(1/7, 7), sides = 1)
  }
  ma
}


############################################################
# 4) FIRST & LAST crossing days for X(t) >= threshold
############################################################
get_cross_days_first_last <- function(X_mat, threshold_y) {
  
  first_day <- apply(X_mat, 2, function(path) {
    idx <- which(path >= threshold_y)
    if (length(idx) == 0) return(NA_integer_)
    idx[1] - 1
  })
  
  last_day <- apply(X_mat, 2, function(path) {
    idx <- which(path >= threshold_y)
    if (length(idx) == 0) return(NA_integer_)
    idx[length(idx)] - 1
  })
  
  list(
    first_day = first_day,
    last_day  = last_day,
    min_first = if (all(is.na(first_day))) NA_integer_ else min(first_day, na.rm = TRUE),
    max_first = if (all(is.na(first_day))) NA_integer_ else max(first_day, na.rm = TRUE),
    min_last  = if (all(is.na(last_day)))  NA_integer_ else min(last_day,  na.rm = TRUE),
    max_last  = if (all(is.na(last_day)))  NA_integer_ else max(last_day,  na.rm = TRUE)
  )
}


############################################################
# 5) VISIBILITY: pick alpha / line width by number of paths
############################################################
pick_band_style <- function(m_draw) {
  alpha <- if (m_draw <= 20) 0.9 else if (m_draw <= 100) 0.4 else if (m_draw <= 500) 0.15 else 0.09
  # alpha = 1
  list(
    col = rgb(0, 0, 0, alpha = alpha),
    lwd = if (m_draw <= 50) 2 else 1
  )
}


############################################################
# 6) GENERIC PLOTTER for moving-average bands
############################################################
plot_ma_band <- function(maX, days, theta, cfg,
                         x_min, x_max, y_min, y_max,
                         ylab, main_title) {
  
  max_paths_to_draw <- cfg$max_paths_to_draw
  
  M <- ncol(maX)
  m_draw <- min(M, max_paths_to_draw)
  cols <- 1:m_draw
  style <- pick_band_style(m_draw)
  
  # Main plot
  matplot(
    x = days, y = maX[, cols, drop = FALSE],
    type = "l", lty = 1, col = style$col, lwd = style$lwd,
    xlab = "Days",
    ylab = ylab,
    main = main_title,
    xlim = c(x_min, x_max),
    ylim = c(y_min, y_max)
  )
  
  # Threshold line + label
  abline(h = theta, lwd = 2, lty = 3, col = "darkred")
  text(
    x = x_min + 0.65 * (x_max - x_min),
    y = theta,
    labels = paste0("threshold = ", format(theta, big.mark = ",")),
    col = "darkred",
    pos = 3
  )
  
  # Crossings
  res <- get_cross_days_first_last(maX, theta)
  
  # First crossings (dark green)
  if (is.finite(res$min_first)) {
    abline(v = res$min_first, lwd = 2, lty = 2, col = "darkgreen")
    text(res$min_first, y_max * 0.95, paste0("min first = ", res$min_first), col = "darkgreen", pos = 4)
  }
  if (is.finite(res$max_first)) {
    abline(v = res$max_first, lwd = 2, lty = 2, col = "darkgreen")
    text(res$max_first, y_max * 0.90, paste0("max first = ", res$max_first), col = "darkgreen", pos = 4)
  }
  
  # Last crossings (purple)
  if (is.finite(res$min_last)) {
    abline(v = res$min_last, lwd = 2, lty = 4, col = "purple")
    text(res$min_last, y_max * 0.85, paste0("min last = ", res$min_last), col = "purple", pos = 2)
  }
  if (is.finite(res$max_last)) {
    abline(v = res$max_last, lwd = 2, lty = 4, col = "purple")
    text(res$max_last, y_max * 0.80, paste0("max last = ", res$max_last), col = "purple", pos = 2)
  }
  
  invisible(res)
}


############################################################
# 7) RUN + COMPUTE maI, maD + PLOT BOTH
############################################################
sim <- simulate_alg1_ID(CFG)

# Daily new
dI <- calc_daily_new(sim$I)
dD <- calc_daily_new(sim$D)

# 7-day moving averages
maI <- calc_ma7(dI)
maD <- calc_ma7(dD)

# Auto y_max (ensure both data + threshold fit)
if (is.na(CFG$y_max_I)) CFG$y_max_I <- max(max(maI, na.rm = TRUE), CFG$theta_I) * 1.05
if (is.na(CFG$y_max_D)) CFG$y_max_D <- max(max(maD, na.rm = TRUE), CFG$theta_D) * 1.05

# Side-by-side plots
par(mfrow = c(1, 2))

resI <- plot_ma_band(
  maX = maI, days = sim$days, theta = CFG$theta_I, cfg = CFG,
  x_min = CFG$x_min_I, x_max = CFG$x_max_I, y_min = CFG$y_min_I, y_max = CFG$y_max_I,
  ylab = "maI(t): 7-day moving avg daily new infections",
  main_title = "Algorithm 1: maI(t)"
)

resD <- plot_ma_band(
  maX = maD, days = sim$days, theta = CFG$theta_D, cfg = CFG,
  x_min = CFG$x_min_D, x_max = CFG$x_max_D, y_min = CFG$y_min_D, y_max = CFG$y_max_D,
  ylab = "maD(t): 7-day moving avg daily new deaths",
  main_title = "Algorithm 1: maD(t)"
)

par(mfrow = c(1, 1))  # reset layout

# Print summaries
cat("maI threshold =", CFG$theta_I, "\n")
cat("  first crossing: min =", resI$min_first, ", max =", resI$max_first, "\n")
cat("  last  crossing: min =", resI$min_last,  ", max =", resI$max_last,  "\n\n")

cat("maD threshold =", CFG$theta_D, "\n")
cat("  first crossing: min =", resD$min_first, ", max =", resD$max_first, "\n")
cat("  last  crossing: min =", resD$min_last,  ", max =", resD$max_last,  "\n")