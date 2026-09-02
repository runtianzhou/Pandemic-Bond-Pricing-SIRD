############################################################
# CUMULATIVE I(t) & D(t) — Algorithm 1 settings (one file)
# - simulate BOTH I and D with correlation rho
# - plot cumulative I and cumulative D (separate plots)
# - horizontal threshold line + min/max crossing vertical lines
# - ONE-PLACE CONFIG at the top (change CFG$M once)
############################################################

#############################
# 0) ONE-PLACE CONFIG
#############################
CFG <- list(
  # -------------------------
  # Simulation parameters
  # -------------------------
  N      = 1104,
  dt     = 1,
  M      = 1,       # <--- change this ONE number to change simulation count everywhere
  I0     = 1,
  D0     = 1,
  KI     = 833e6,
  KD     = 0.4e6,
  gI     = 0.320,
  gD     = 0.320,
  sigmaI = 0.1,
  sigmaD = 0.1,
  rho    = 0.5,
  seed   = NULL,       # set 42 for reproducible results
  
  # -------------------------
  # Threshold lines (reference levels on cumulative plots)
  # -------------------------
  threshold_I = 5000,
  threshold_D = 2500,
  
  # -------------------------
  # Plot zoom controls (I plot)
  # -------------------------
  x_min_I = 0,
  x_max_I = 200,       # set 1104 for full horizon
  y_min_I = 0,
  y_max_I = NA,
  
  # -------------------------
  # Plot zoom controls (D plot)
  # -------------------------
  x_min_D = 0,
  x_max_D = 200,       # set 1104 for full horizon
  y_min_D = 0,
  y_max_D = NA,       # adjust as you like
  
  # -------------------------
  # Plot performance control
  # -------------------------
  max_paths_to_draw = 5000
)


############################################################
# 1) SIMULATE CUMULATIVE I(t), D(t) (Euler–Maruyama + correlation)
############################################################
simulate_alg1_ID <- function(cfg) {
  
  # Unpack
  N <- cfg$N; dt <- cfg$dt; M <- cfg$M
  I0 <- cfg$I0; D0 <- cfg$D0
  KI <- cfg$KI; KD <- cfg$KD
  gI <- cfg$gI; gD <- cfg$gD
  sigmaI <- cfg$sigmaI; sigmaD <- cfg$sigmaD
  rho <- cfg$rho
  seed <- cfg$seed
  
  # Checks
  if (N < 1) stop("N must be >= 1")
  if (dt <= 0) stop("dt must be > 0")
  if (M < 1) stop("M must be >= 1")
  if (KI <= 0 || KD <= 0) stop("KI and KD must be > 0")
  if (rho < -1 || rho > 1) stop("rho must be in [-1, 1]")
  
  set.seed(seed)
  
  # Storage
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
    
    # Correlated increments:
    # dBI = dW1
    # dBD = rho*dW1 + sqrt(1-rho^2)*dW2
    dBI <- dW1
    dBD <- rho * dW1 + sqrt_1mrho2 * dW2
    
    I_prev <- I[t - 1, ]
    D_prev <- D[t - 1, ]
    
    # Stochastic logistic Euler–Maruyama updates
    I_new <- I_prev + I_prev * (1 - I_prev / KI) * (gI * dt + sigmaI * dBI)
    D_new <- D_prev + D_prev * (1 - D_prev / KD) * (gD * dt + sigmaD * dBD)
    
    # Clamp to [0, K] for stability
    I[t, ] <- pmin(pmax(I_new, 0), KI)
    D[t, ] <- pmin(pmax(D_new, 0), KD)
  }
  
  list(I = I, D = D, days = 0:N)
}


############################################################
# 2) INTERSECTION DAYS: first day X(t) >= threshold
############################################################
get_cross_days_minmax <- function(X_mat, threshold_y) {
  
  cross_day <- apply(X_mat, 2, function(path) {
    idx <- which(path >= threshold_y)
    if (length(idx) == 0) return(NA_integer_)
    idx[1] - 1  # row 1 is day 0
  })
  
  list(
    cross_day = cross_day,
    min_day   = if (all(is.na(cross_day))) NA_integer_ else min(cross_day, na.rm = TRUE),
    max_day   = if (all(is.na(cross_day))) NA_integer_ else max(cross_day, na.rm = TRUE)
  )
}


############################################################
# 3) Helper: choose color opacity/line width based on paths drawn
############################################################
pick_band_style <- function(m_draw) {
  alpha <- if (m_draw <= 20) 0.9 else if (m_draw <= 100) 0.4 else if (m_draw <= 500) 0.15 else 0.03
  list(
    col = rgb(0, 0, 0, alpha = alpha),
    lwd = if (m_draw <= 50) 2 else 1
  )
}


############################################################
# 4) Generic plot function for cumulative paths
############################################################
plot_cumulative_with_threshold_and_crossings <- function(
    X, days,
    threshold_y,
    x_min, x_max, y_min, y_max,
    max_paths_to_draw = 5000,
    ylab = "Cumulative",
    main = "",
    threshold_col = "darkred",
    cross_col = "darkgreen"
) {
  
  M <- ncol(X)
  m_draw <- min(M, max_paths_to_draw)
  cols <- 1:m_draw
  
  style <- pick_band_style(m_draw)
  
  # Main band plot
  matplot(
    x = days, y = X[, cols, drop = FALSE],
    type = "l", lty = 1, col = style$col, lwd = style$lwd,
    xlab = "Days", ylab = ylab,
    main = main,
    xlim = c(x_min, x_max),
    ylim = c(y_min, y_max)
  )
  
  # Threshold line + label
  abline(h = threshold_y, lwd = 2, lty = 3, col = threshold_col)
  text(
    x = x_min + 0.65 * (x_max - x_min),
    y = threshold_y,
    labels = paste0("threshold = ", format(threshold_y, big.mark = ",")),
    col = threshold_col,
    pos = 3
  )
  
  # Min/max crossing days across ALL paths
  res <- get_cross_days_minmax(X, threshold_y)
  min_day <- res$min_day
  max_day <- res$max_day
  
  # Vertical lines + labels
  if (is.finite(min_day)) {
    abline(v = min_day, lwd = 2, lty = 2, col = cross_col)
    text(min_day, y_max * 0.95, paste0("min day = ", min_day), col = cross_col, pos = 4)
  }
  if (is.finite(max_day)) {
    abline(v = max_day, lwd = 2, lty = 2, col = cross_col)
    text(max_day, y_max * 0.90, paste0("max day = ", max_day), col = cross_col, pos = 4)
  }
  
  invisible(res)
}


############################################################
# 5) RUN + PLOT BOTH
############################################################
sim <- simulate_alg1_ID(CFG)

# Auto-set Y max if left as NA
if (is.na(CFG$y_max_I)) {
  CFG$y_max_I <- max(sim$I, na.rm = TRUE) * 1
}

if (is.na(CFG$y_max_D)) {
  CFG$y_max_D <- max(sim$D, na.rm = TRUE) * 1
}

# # (sets the plotting layout for the next plots)
par(mfrow = c(1, 2))   # 1 row, 2 columns (side-by-side)


# Plot cumulative I(t)
resI <- plot_cumulative_with_threshold_and_crossings(
  X = sim$I, days = sim$days,
  threshold_y = CFG$threshold_I,
  x_min = CFG$x_min_I, x_max = CFG$x_max_I,
  y_min = CFG$y_min_I, y_max = CFG$y_max_I,
  max_paths_to_draw = CFG$max_paths_to_draw,
  ylab = "Cumulative infected I(t)",
  main = "Algorithm 1: Cumulative I(t)"
)

cat("I(t) threshold =", CFG$threshold_I, "\n")
cat("  min crossing day:", resI$min_day, "\n")
cat("  max crossing day:", resI$max_day, "\n\n")

# Plot cumulative D(t)
resD <- plot_cumulative_with_threshold_and_crossings(
  X = sim$D, days = sim$days,
  threshold_y = CFG$threshold_D,
  x_min = CFG$x_min_D, x_max = CFG$x_max_D,
  y_min = CFG$y_min_D, y_max = CFG$y_max_D,
  max_paths_to_draw = CFG$max_paths_to_draw,
  ylab = "Cumulative deaths D(t)",
  main = "Algorithm 1: Cumulative D(t)"
)

cat("D(t) threshold =", CFG$threshold_D, "\n")
cat("  min crossing day:", resD$min_day, "\n")
cat("  max crossing day:", resD$max_day, "\n")